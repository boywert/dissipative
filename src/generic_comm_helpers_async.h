#ifdef GENERIC_ASYNC

#include <stdio.h>
#include <stdlib.h>

#ifdef DISABLE_MEMORY_MANAGER
#error "GENERIC_ASYNC cannot be used together with DISABLE_MEMORY_MANAGER"
#endif

#ifdef USE_SUBCOMM_COMMUNICATOR
#define MYCOMMUNICATOR SubComm
#define MyThisTask  SubThisTask
#define MyNTask     SubNTask
#define MyTagOffset SubTagOffset
#else
#define MYCOMMUNICATOR MPI_COMM_WORLD
#define MyThisTask  ThisTask
#define MyNTask     NTask
#define MyTagOffset TagOffset
#endif

#define EXTRA_SPACE 16384

#define POLLINGINTERVAL  255

#define MIN_PRIMARY_FRACTION_BEFORE_SERVING_OF_REQUESTS   0.2

#define MEMORY_FRACTION   0.333
#define EXTRAFAC          3                             /* allow this many more work requests than processor number to be concurrently processed */


typedef struct datanodelist datanodelist;
typedef struct data_partlist data_partlist;


static size_t ExportSpace;
static size_t MinSpace;
static int NextParticle;
static int Nexport, Nimport;
static int NexportNodes, NimportNodes;
static long long SumNexport;
static int *NodeDataIn;
static int *NodeDataGet;

static char callorigin[1000];

static MPI_Request *result_recv_requests;
static MPI_Request *result_send_requests;
static MPI_Request *requests;

static int n_result_recv_requests;
static int n_result_send_requests;
static int n_requests;
static int n_send_buffers;

static int *compl_list;

static data_out **dataresult_list;



#define generic_set_MaxNexport(...) { generic_set_info(__FUNCTION__, __FILE__, __LINE__); }

/* this function determines how much buffer space we may use based on the memory that is locally still free,
 * and it computes how much memory may at most be needed to process a single particle. We will only continue with a particle
 * if this can still be safely processed.
 */
static void generic_set_info(const char *func, const char *file, int line)
{
  ExportSpace = MEMORY_FRACTION * FreeBytes;  /* we just grab at most this fraction of the still available memory here */

  ExportSpace /= NUM_THREADS;
  ExportSpace -= NumPart * sizeof(int);  /* to account for the neighbor list buffer that every thread allocated */

  /* make the size a multiple both of data_partlist and datanodelist */
  ExportSpace /= (sizeof(data_partlist) * sizeof(datanodelist));
  ExportSpace *= (sizeof(data_partlist) * sizeof(datanodelist));

  MinSpace = (MyNTask - 1) * (sizeof(data_partlist) + sizeof(data_in) + sizeof(data_out)) + NTopleaves * (sizeof(datanodelist) + sizeof(int));

  sprintf(callorigin, "%s|%d|", file, line);

#ifdef VERBOSE
  mpi_printf("GENERIC-ASYNC: function %s(), file %s, line %d: MinSpace = %g MB  NTopleaves = %d  ExportSpace = %g MB\n", func, file, line,
             MinSpace * TO_MBYTE_FAC, NTopleaves, ExportSpace * TO_MBYTE_FAC);
#endif

  if(ExportSpace < MinSpace)
    terminate("Bummer. Can't even safely process a single particle for the available memory. FreeBytes=%lld  ExportSpace=%lld  MinSpace=%lld  MyNTask=%d  NTopleaves=%d",
              (long long) FreeBytes, (long long) ExportSpace, (long long) MinSpace, MyNTask, NTopleaves);
}


/* this function does the memory allocation at the beginning of a loop over the remaining local particles.
 * The fields PartList[] and NodeList[] share the buffer space of size "ExportSpace" (in bytes).
 * Here PartList will be filled in from the beginning, while NodeList will be filled in from the end.
 * Since we do not know a priory the relative share of these two fields, we can make optimum use of
 * the available space in this way.
 */
static void generic_alloc_partlist_nodelist_ngblist_threadbufs(void)
{
  Thread[0].Ngblist = (int *) mymalloc_g("Ngblist", NumPart * (NUM_THREADS * sizeof(int)));
  for(int i = 1; i < NUM_THREADS; i++)
    Thread[i].Ngblist = Thread[i - 1].Ngblist + NumPart;

  Thread[0].R2list = (double *) mymalloc_g("R2list", NumPart * (NUM_THREADS * sizeof(double)));
  for(int i = 1; i < NUM_THREADS; i++)
    Thread[i].R2list = Thread[i - 1].R2list + NumPart;
 }

/* the corresponding deallocation routine
 */
static void generic_free_partlist_nodelist_ngblist_threadbufs(void)
{
  myfree(Thread[0].R2list);
  myfree(Thread[0].Ngblist);

  for(int i = NUM_THREADS - 1; i >= 0; i--)
    {
      Thread[i].R2list = NULL;
      Thread[i].Ngblist = NULL;
    }

}

static void generic_initialize_counters(void)
{
  for(int i = 0; i < NUM_THREADS; i++)
    {
      Thread[i].Nexport = 0;
      Thread[i].NexportNodes = 0;
      Thread[i].ExportSpace = ExportSpace;
      Thread[i].InitialSpace = ExportSpace;
      Thread[i].ItemSize = (sizeof(data_partlist) + sizeof(data_in) + sizeof(data_out));

      if(i > 0)
	Thread[i].PartList = (struct data_partlist *) ((char *) Thread[i - 1].PartList + ExportSpace);

      Thread[i].Exportflag = Exportflag + i * ((((MyNTask - 1) / 16) + 1) * 16);
    }
}


static size_t generic_consolidate_partlist_storage(void)
{
  size_t sumgaps = 0;

  for(int i = 0; i < NUM_THREADS; i++)
    {
      memmove((char *) Thread[i].PartListOld - sumgaps, (char *) Thread[i].PartList, Thread[i].Nexport * sizeof(data_partlist));

      size_t gap = ExportSpace - Thread[i].Nexport * sizeof(data_partlist) -  Thread[i].NexportNodes * sizeof(datanodelist);

      /* make the gap a multiple both of data_partlist and datanodelist */
      gap /= (sizeof(data_partlist) * sizeof(datanodelist));
      gap *= (sizeof(data_partlist) * sizeof(datanodelist));

      memmove((char *) Thread[i].PartListOld + Thread[i].InitialSpace - Thread[i].NexportNodes * sizeof(datanodelist) - sumgaps - gap,
              (char *) Thread[i].PartList + Thread[i].InitialSpace - Thread[i].NexportNodes * sizeof(datanodelist),
              Thread[i].NexportNodes * sizeof(datanodelist));

      Thread[i].PartListOld = (struct data_partlist *) ((char *) Thread[i].PartListOld - sumgaps);
      Thread[i].InitialSpace -= gap;

      sumgaps += gap;
    }

  return NUM_THREADS * ExportSpace - sumgaps;
}


static void generic_prepare_export_counts(void)
{
  for(int j = 0; j < MyNTask; j++)
    {
      Send[j].Count = 0;
      Send[j].CountNodes = 0;
    }

  Nexport = 0;
  NexportNodes = 0;

  for(int i = 0; i < NUM_THREADS; i++)
    {
      for(int j = 0; j < Thread[i].Nexport; j++)
        Send[Thread[i].PartList[j].Task].Count++;

      datanodelist *nodelist = (datanodelist *) (((char *) Thread[i].PartList) + Thread[i].InitialSpace);

      for(int j = 0; j < Thread[i].NexportNodes; j++)
        Send[nodelist[-1 - j].Task].CountNodes++;

      Nexport += Thread[i].Nexport;
      NexportNodes += Thread[i].NexportNodes;
    }

  SumNexport += Nexport;
}



/* initialize offset tables that we need for the communication
 */
static void generic_prepare_export_offsets(void)
{
  Send_offset[0] = 0;
  Send_offset_nodes[0] = 0;

  for(int j = 1; j < MyNTask; j++)
    {
      Send_offset[j] = Send_offset[j - 1] + Send[j - 1].Count;
      Send_offset_nodes[j] = Send_offset_nodes[j - 1] + Send[j - 1].CountNodes;
    }
}


/* organize the particle and node data for export in contiguous memory regions
 */
static void generic_prepare_particle_data_for_export(void)
{
  int *rel_node_index = (int *) mymalloc_g("rel_node_index", MyNTask * sizeof(int));

  for(int j = 0; j < MyNTask; j++)
    {
      Send[j].Count = 0;
      Send[j].CountNodes = 0;
      rel_node_index[j] = 0;
    }

  for(int i = 0; i < NUM_THREADS; i++)
    {
      datanodelist *nodelist = (datanodelist *) (((char *) Thread[i].PartList) + Thread[i].InitialSpace);

      for(int j = 0, jj = 0; j < Thread[i].Nexport; j++)
        {
          int task = Thread[i].PartList[j].Task;
          int off = Send_offset[task] + Send[task].Count++;

          int target = Thread[i].PartList[j].Index;

          particle2in(&DataIn[off], target, rel_node_index[task]);

          if(j < Thread[i].Nexport - 1)
            if(Thread[i].PartList[j].Index == Thread[i].PartList[j + 1].Index)
              continue;

          while(jj < Thread[i].NexportNodes && Thread[i].PartList[j].Index == nodelist[-1 - jj].Index)
            {
              task = nodelist[-1 - jj].Task;
              off = Send_offset_nodes[task] + Send[task].CountNodes++;

              NodeDataIn[off] = nodelist[-1 - jj].Node;

              rel_node_index[task]++;
              jj++;
            }
        }
    }

  myfree(rel_node_index);
}

/* driver routine to process the results that we obtained for a particle from a remote processor
 * by working on it with the supplied out2particle() routine
 */
static void generic_add_results_to_local(void)
{
  for(int j = 0; j < MyNTask; j++)
    Send[j].Count = 0;

  for(int i = 0; i < NUM_THREADS; i++)
    for(int j = 0; j < Thread[i].NexportOld; j++)
      {
        int task = Thread[i].PartListOld[j].Task;
        int off = Send_offset[task] + Send[task].Count++;

        int target = Thread[i].PartListOld[j].Index;

        out2particle(&DataOut[off], target, MODE_IMPORTED_PARTICLES);
      }
}

#ifndef OMIT_GENERIC_GET_NUMNODES

/* this function is called in the actual tree walk routine to find out how the number and
 * starting index of the section in the node-list that needs to be processed for the imported particle
 */
static void generic_get_numnodes(int target, int *numnodes, int **firstnode)
{
  if(target == Nimport - 1)
    *numnodes = NimportNodes - DataGet[target].Firstnode;
  else
    *numnodes = DataGet[target + 1].Firstnode - DataGet[target].Firstnode;

  *firstnode = &NodeDataGet[DataGet[target].Firstnode];
}

#endif


/* calculate how many space we need to allocate to safely process a certain number of
 * nodes and particles that are imported.
 */
static size_t generic_calc_import_storage(int nimport, int nimportnodes)
{
  size_t needed = nimport * sizeof(data_in) + nimportnodes * sizeof(int) + nimport * sizeof(data_out);

  /* add some extra space to not go to the last byte */
  needed += EXTRA_SPACE;

  return needed;
}

static int generic_sort_ints(const void *p1, const void *p2)
{
  if(*(int *)p1 < *(int *)p2)
    return -1;
  else if(*(int *)p1 > *(int *)p2)
    return +1;

  return 0;
}

/* This function removes 'count' entries from requests[] of length 'n', with indices given in list[].
 * Note, the array list[] must be of length n as well as we are using the extra storage this gives */


static void generic_remove_from_request_list(int n, MPI_Request *requests, int count, int *list)
{
  if(count == n)
    return;

  if(count > 1)
    qsort(list, count, sizeof(int), generic_sort_ints);
 
  int *newplace = (int *) mymalloc_g("newplace", sizeof(int) * n);

  for(int i = 0; i < count; i++)
    {
      int idx = list[i];
      if(idx >= n)
        idx = newplace[idx];

      requests[idx] = requests[n - 1];
      newplace[n - 1] = idx;
      n--;
    }

  myfree(newplace);
}

static void generic_polling(void)
{
  int ncompleted;

  if(n_result_recv_requests)
    {
      MPI_Testsome(n_result_recv_requests, result_recv_requests, &ncompleted, compl_list, MPI_STATUSES_IGNORE);

      if(ncompleted)
        {
          generic_remove_from_request_list(n_result_recv_requests, result_recv_requests, ncompleted, compl_list);
          n_result_recv_requests -= ncompleted;
        }
    }

  if(n_requests)
    {
      MPI_Testsome(n_requests, requests, &ncompleted, compl_list, MPI_STATUSES_IGNORE);

      if(ncompleted)
        {
          generic_remove_from_request_list(n_requests, requests, ncompleted, compl_list);
          n_requests -= ncompleted;
        }
    }

  if(n_result_send_requests)
    {
      MPI_Testsome(n_result_send_requests, result_send_requests, &ncompleted, compl_list, MPI_STATUSES_IGNORE);

      if(ncompleted)
        {
          generic_remove_from_request_list(n_result_send_requests, result_send_requests, ncompleted, compl_list);
          n_result_send_requests -= ncompleted;
        }
    }
}


static int generic_polling_primary(int count, int Nforces)
{
  if(count > MIN_PRIMARY_FRACTION_BEFORE_SERVING_OF_REQUESTS * Nforces + 1)
    {
      int flag;

      /* check whether we got a work request */
      MPI_Iprobe(MPI_ANY_SOURCE, TAG_N + MyTagOffset, MYCOMMUNICATOR, &flag, MPI_STATUS_IGNORE);

      if(flag)
        return flag;
    }

  generic_polling();

  return 0;
}


static void generic_polling_secondary(void)
{
  generic_polling();
}



static int generic_check_for_work_request(void (*kernel_imp) (void))
{
  int flag = 0;
  MPI_Status status;

  /* check whether we got a work request */
  
  if(n_send_buffers < EXTRAFAC * MyNTask)
    MPI_Iprobe(MPI_ANY_SOURCE, TAG_N + MyTagOffset, MYCOMMUNICATOR, &flag, &status);

  if(flag)
    {
      int source = status.MPI_SOURCE;
      int tag = status.MPI_TAG;
      
      MPI_Recv(&Recv[source], sizeof(struct send_recv_counts), MPI_BYTE, source, tag, MYCOMMUNICATOR, MPI_STATUS_IGNORE);
      
      Nimport = Recv[source].Count;
      NimportNodes = Recv[source].CountNodes;
      
      if(generic_calc_import_storage(Nimport, NimportNodes) > FreeBytes)
        terminate("Sorry, we have a serious problem: %g MB byte are free byut we need %g for this import",
            FreeBytes * TO_MBYTE_FAC, generic_calc_import_storage(Nimport, NimportNodes) * TO_MBYTE_FAC);

      /* now allocated the import and results buffers */
      DataResult = (data_out *) mymalloc_g("DataResult", Nimport * sizeof(data_out));
      DataGet = (data_in *) mymalloc_g("DataGet", Nimport * sizeof(data_in));
      NodeDataGet = (int *) mymalloc_g("NodeDataGet", NimportNodes * sizeof(int));
      
      /* get the particles */
      MPI_Recv(DataGet, Nimport * sizeof(data_in), MPI_BYTE, source, TAG_PART_DATA + MyTagOffset, MYCOMMUNICATOR, MPI_STATUS_IGNORE);
      
      /* get the nodes */
      MPI_Recv(NodeDataGet, NimportNodes, MPI_INT, source, TAG_NODE_DATA + MyTagOffset, MYCOMMUNICATOR, MPI_STATUS_IGNORE);
      
      /* now do the actual work for the imported points */
      kernel_imp();
      
      myfree(NodeDataGet);
      myfree(DataGet);
      
      /* send the results */

      dataresult_list[n_send_buffers++] = DataResult;
      MPI_Issend(DataResult, Recv[source].Count * sizeof(data_out), MPI_BYTE, source, TAG_RESULTS + MyTagOffset, MYCOMMUNICATOR, &result_send_requests[n_result_send_requests++]);
    }
  return flag;
}



/* Implements a repeated loop over the local particles in the list, processing them with the local kernel function,
 * until we're done or the export buffer is full. Then we exchange the data, and process the imported ones with the provided kernel.
 * We repeat if neeed until all processors are done.
 */
static int generic_comm_pattern(int nactive, void (*kernel_loc) (void), void (*kernel_imp) (void))
{
  int iter = 0;

  SumNexport = 0;               /* can be queried as a book-keeping variable */

  NextParticle = 0;             /* first particle index for this task */

  MyTagOffset++;
  if(MyTagOffset > 30000)  /* the MPI Standard defines 32767 as the minimum value for the largest tag that must be allowed */
    MyTagOffset = 0;

  /* variables for non-blacking barrier */
  int nLevels = my_fls(MyNTask - 1);
  int received_levels = 0, sent_levels = 0;
  int *stagelist = (int *) mymalloc_g("stagelist", nLevels * sizeof(int) );
  for(int j = 0; j < nLevels; j++)
    stagelist[j] = j;


  MPI_Request *level_requests = (MPI_Request *) mymalloc_g("level_requests", nLevels * sizeof(MPI_Request));

  compl_list = (int *) mymalloc_g("compl_list", imax(3, EXTRAFAC) * MyNTask * sizeof(int));

  /* queue A: collects the sending of work requests */
  requests = (MPI_Request *) mymalloc_g("requests", 3 * MyNTask * sizeof(MPI_Request));
  n_requests = 0;

  /* queue B: collects the receiving of partial results */
  result_recv_requests = (MPI_Request *) mymalloc_g("result_recv_requests", MyNTask * sizeof(MPI_Request));
  n_result_recv_requests = 0;

  /* queue C: collects the sending of partial results */
  result_send_requests = (MPI_Request *) mymalloc_g("result_send_requests", EXTRAFAC * MyNTask * sizeof(MPI_Request));
  n_result_send_requests = 0;
  n_send_buffers = 0;

  dataresult_list = (data_out **) mymalloc_g("dataresult_list", EXTRAFAC * MyNTask * sizeof(data_out *));



  /* allocate buffers to arrange communication */
  generic_alloc_partlist_nodelist_ngblist_threadbufs();

  for(int i = 0; i < NUM_THREADS; i++)
    {
      Thread[i].PartList = NULL;
      Thread[i].Nexport = 0;
    }

  int barrier_active = 0;

  while(1)
    {
      /* let's temporarily store PartList/Nexport in PartListOld/NexportOld to be able to delay the adding in
         of partial results, so that we can anticipate a call of the local kernel function to already start
         processing local particles */
      for(int i = 0; i < NUM_THREADS; i++)
        {
          Thread[i].PartListOld =  Thread[i].PartList;
          Thread[i].NexportOld = Thread[i].Nexport;
        }

      Thread[0].PartList = (struct data_partlist *) mymalloc_g("PartList", NUM_THREADS * ExportSpace);

      generic_initialize_counters();

      /* do local particles if we are not done yet */
      if(NextParticle < nactive)
        {
          kernel_loc();
          iter++;
        }

      /* if some results we sent in previous round haven't been picked up yet, need to wait for them now */
      while(n_result_send_requests || n_result_recv_requests || n_requests)
        {
          generic_check_for_work_request(kernel_imp);

          generic_polling();
        }

      /* add the results to the local particles (use OldPartList for this) */
      generic_add_results_to_local();

      /* set up (new) Sendcount table */
      generic_prepare_export_counts();

      /* prepare offsets in export tables */
      generic_prepare_export_offsets();

      if(Thread[0].PartListOld)
	{
          /* We do a really dirty trick now in order to avoid having to move full memory buffers. We exploit that
           * we know that DataIn, NodeDataIn, etc. are allocated at this point right after PartListOld, and that
           * none of this data behind this block is needed any more.
           * So even if the storage needed for PartList is larger than the space originally at PartListOld,
           * nothing bad happends. This method is OK for our stack-wise allocation, but wouldn't work for the ordinary
           * independent allocation from the system heap, this is why we need to disallow DISABLE_MEMORY_MANAGER here.
           */

	  /* reinflate the other thread pointers so that PartListOld looks like freshly allocated */
	  for(int i = 1; i < NUM_THREADS; i++)
	    Thread[i].PartListOld = (struct data_partlist *) ((char *) Thread[i - 1].PartListOld + ExportSpace);

	  /* now arrange the data in contiguous blocks in PartListOld (this also moves it forward on the stack) */
	  size_t new_size = generic_consolidate_partlist_storage();

	  /* free all the send buffer blocks, plus one more, which is PartList. This is here done in this way because
	   * some of the send data blocks may appear before, some after the allocation of PartList
	   */
	  for(int i = 0; i <= n_send_buffers; i++)
	    myfree(myfree_query_last_block());

	  /* now call our PartListOld pointers PartList again */
	  for(int i = 0; i < NUM_THREADS; i++)
	    Thread[i].PartList = Thread[i].PartListOld;

	  /* explicitly mark the data block buffers as unallocated (just for safety so that we are guaranteed not to use them by accident) */
	  for(int i = n_send_buffers - 1; i >= 0; i--)
	    dataresult_list[i] = NULL;

	  n_send_buffers = 0;

	  myfree(DataOut);    
	  myfree(NodeDataIn); 
	  myfree(DataIn);     

	  /* finally, we can shrink PartList, which is now the final block, to the new size */
	  myrealloc(Thread[0].PartList, new_size);
	}
      else
	{
          /* here we treat the first iteration when no result blocks are present yet */
	  for(int i = 0; i < NUM_THREADS; i++)
	    Thread[i].PartListOld = Thread[i].PartList;

	  /* now arrange the data in contiguous blocks  */
	  size_t new_size = generic_consolidate_partlist_storage();

	  myrealloc(Thread[0].PartList, new_size);
	}


      /* now prepare data from this round*/

      /* allocate particle data buffers */
      DataIn = (data_in *) mymalloc_g("DataIn", Nexport * sizeof(data_in));
      NodeDataIn = (int *) mymalloc_g("NodeDataIn", NexportNodes * sizeof(int));
      DataOut = (data_out *) mymalloc_g("DataOut", Nexport * sizeof(data_out));

      /* prepare particle data for export */
      generic_prepare_particle_data_for_export();


      /* post our export work requests, as well as our result-receive requests */
      for(int ngrp = 0; ngrp < (1 << PTask); ngrp++)
        {
          int j = MyThisTask ^ ngrp;

          Recv[j].Count = 0;
          Recv[j].CountNodes = 0;

          if(j < MyNTask)
            if(Send[j].Count > 0)
              {
                MPI_Issend(&Send[j], sizeof(struct send_recv_counts), MPI_BYTE, j, TAG_N + MyTagOffset, MYCOMMUNICATOR, &requests[n_requests++]);
                MPI_Issend(&DataIn[Send_offset[j]], Send[j].Count * sizeof(data_in), MPI_BYTE, j, TAG_PART_DATA + MyTagOffset, MYCOMMUNICATOR, &requests[n_requests++]);
                MPI_Issend(&NodeDataIn[Send_offset_nodes[j]], Send[j].CountNodes, MPI_INT, j, TAG_NODE_DATA + MyTagOffset, MYCOMMUNICATOR, &requests[n_requests++]);
                MPI_Irecv(&DataOut[Send_offset[j]], Send[j].Count * sizeof(data_out), MPI_BYTE, j, TAG_RESULTS + MyTagOffset, MYCOMMUNICATOR, &result_recv_requests[n_result_recv_requests++]);
            }
        }

      while(1)  /* polling loop */
        {
          /* first, check whether we got a work request, otherwise look for barrier signals */
          if(!generic_check_for_work_request(kernel_imp))
            {
	      int flag = 0;
	      MPI_Status status;

              /* now check whether we have one of the barrier markers  */
              MPI_Iprobe(MPI_ANY_SOURCE, TAG_BARRIER + MyTagOffset, MYCOMMUNICATOR, &flag, &status);

              if(flag)
                {
                  int source = status.MPI_SOURCE;
                  int tag = status.MPI_TAG;

                  int stage;
                  MPI_Recv(&stage, 1, MPI_INT, source, tag, MYCOMMUNICATOR, MPI_STATUS_IGNORE);
                  received_levels |= (1 << stage);
                }
              else
                {
		  /* looks like we are done for this task, hence we may start the asynchronous barrier */
                  if(NextParticle >= nactive && n_result_recv_requests == 0)
                    barrier_active = 1;

                  if(barrier_active == 1)
                    {
                      for(int stage = 0; stage < nLevels; stage++)
                        if(!(sent_levels & (1 << stage)))
                          {
                            int mask = ((1 << stage) - 1);

                            if((mask & received_levels) == mask)
                              {
                                sent_levels |= (1 << stage);

                                int target = (MyThisTask + (1 << stage)) % MyNTask;

                                MPI_Issend(&stagelist[stage], 1, MPI_INT, target, TAG_BARRIER + MyTagOffset, MYCOMMUNICATOR, &level_requests[stage]);
                              }
                          }
                    }
                  break;
                }
            }
        }

      if(received_levels == ((1 << nLevels) - 1) && sent_levels == ((1 << nLevels) - 1))
	{
	  /* make sure that all our barrier sends have really been dealt with and have been picked up (note that we are going to free stagelist) */
	  int flag;
	  MPI_Testall(nLevels, level_requests, &flag, MPI_STATUSES_IGNORE);
	  if(flag)
	    break;
	}
    }

  /* make sure that any open sends of partial results have really been picked up */
  if(n_result_send_requests)
    MPI_Waitall(n_result_send_requests, result_send_requests, MPI_STATUSES_IGNORE);

  /* free any remaining allocated send buffers */
  for(int i = n_send_buffers - 1; i >= 0; i--)
    myfree(dataresult_list[i]);

  n_send_buffers = 0;

  myfree(DataOut);
  myfree(NodeDataIn);
  myfree(DataIn);
  myfree(Thread[0].PartList);
  
  /* free the rest of the buffers */
  generic_free_partlist_nodelist_ngblist_threadbufs();

  myfree(dataresult_list);
  myfree(result_send_requests);
  myfree(result_recv_requests);
  myfree(requests);
  myfree(compl_list);
  myfree(level_requests);
  myfree(stagelist);

  return iter;
}

#endif

