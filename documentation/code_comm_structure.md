
Generic Communication Pattern
=============================

Here we provide some detailed information on how to use the new "generic" communication helpers in the 
parallel tree walk that were added to Arepo in r30950. *This description originally from the post 
titled "new comm structure introduced in r30950" from Volker.*

Recall first the basic way how parallel tree walks (most often for working on a local neighbourhood) 
algorithmically work in Arepo. For a given set of target particles (or cells/coordinates), first a 
local ("primary") tree walk is done to find all interacting particles stored on the local MPI rank. 
At the same time, this local tree walk notices which other domains (and hence MPI ranks) may have 
further interacting particles/cells, and at which node of the tree on the remote processor one needs 
to continue the tree walk. This information is written into an export buffer. After the primary tree 
walks have finished (either because all local particles are done or the export buffer is full), the 
content of the export buffers is communicated to the right target processors and a corresponding 
import buffer is filled on each MPI rank. For these imported points, a "secondary" tree walk is now 
started, just for those branches of the local tree that could not be accessed by the remote processor. 
The additional interactions found are processed accordingly, and finally, the result of this is sent 
back to the requesting MPI rank and is added to (or combined with) the result of the local tree walk. 
The whole procedure may have to be repeated until all particles are processed if the export buffer 
full condition had been encountered.

This pattern is followed in several dozen routines in Arepo, and each time a lot of very similar code 
is executed (and was repeated in the source) that deals with the bookkeeping involved. This makes it 
not only quite error prone but also obscures what's actually done in a particular routine, because 
sometimes the bookkeeping code could in fact be longer than the non-trivial actions of the routine. 
Towards the possibility of full thread support in Arepo, the old communication pattern was changed in 
r30950, at which point things were simplified by moving most of the bookkeeping code to a header file. 
This is then in part used as a poor man's templating mechanism. The way this works is as follows:


**[1]** With this mechanism, each \*.c-source file can implement only one parallel tree walk. At the 
beginning of the file, two static structure types with two variables each have to be defined::

    static struct data_in
    {
      /* your fields go here */
      
      int Firstnode;
    }
     *DataIn, *DataGet;

    static struct data_out
    {
      /* your fields go here */

    } *DataResult, *DataOut;

The names need to be kept exactly as given here. In the ``data_in`` structure you should add all input 
fields that are needed to carry out the relevant part of the tree walk on a remote processor, while in 
the ``data_out`` structure you should add fields required to store the corresponding result. For example, 
for a SPH-like density calculation, this would be the position and smoothing length in ``data_in``, and 
the (partial) density in ``data_out``.


**[2]** Next, you need to define two static functions at the beginning of the file::

    static void particle2in(struct data_in *in, int i, int firstnode)
    {
      /* your field copy statements go here */

      in->Firstnode = firstnode;
    }

    static void out2particle(struct data_out *out, int i, int mode)
    {
      if(mode == MODE_LOCAL_PARTICLES)      /* initial store */
        {
          /* your statement to store a local result go here */
        }
      else                                  /* combine */
        {
          /* your statements to add a remote result to a local one go here */
        }
    }

Again, this template needs to be kept unchanged, in particular, the statement ``in->Firstnode = firstnode;`` 
always needs to there. The purpose of the first function is to fill in the data that needs to be sent to a 
foreign processor into the structure pointed to by ``in``, while the purpose of the second function is to 
deal with the results stored in the structure pointed to by ``out``. In both cases, ``in`` and ``out`` are 
associated with a particle indexed by the argument ``i``. For the example of an SPH-like density 
computation, the body of the first function may, e.g., look like::

   in->Pos[0] = P[i].Pos[0];
   in->Pos[1] = P[i].Pos[1];
   in->Pos[2] = P[i].Pos[2];
   in->Hsml = SphP[i].Hsml;
   in->Firstnode = firstnode;

while the body of the second function could look like this::

   if(mode == MODE_LOCAL_PARTICLES)   
     {
       SphP[i].Density = out->Density;
     }
   else                               
     {
       SphP[i].Density += out->Density;
     }


**[3]** *After* these preliminaries, you add::

    #include "generic_comm_helpers.h"

So this header file is included not at the beginning of the file, but after the above definitions 
and code snippets.


**[4]** Let's now discuss the code for the actual parallel tree walk computation in this framework. 
We'll assume that this is best split into a main function (which loops over all local target 
particles), and one that does a tree walk, either for a local particle, or for an imported one. The 
second one is a purely local function that should be declared static in the file, let's call that 
``tree_evaluate`` here. We then add at the beginning of the file::

    static void tree_evaluate(int target, int mode, int threadid);

Our main function will be called externally, and could for example be one that is supposed to inject 
some feedback from stars in their local neighbourhood. We might call this::

    void inject_feedback_from_stars(void);

and add a corresponding prototype in ``proto.h``.

The generic implementation of the main function should be similar to this::

    void inject_feedback_from_stars(void)
    {
       generic_set_MaxNexport();

       NextParticle = 0;   /* first particle index for this task */

       do
          {
              generic_alloc_partlist_nodelist_ngblist_threadbufs();

    #pragma omp parallel private(i, idx)
              {
                int j, threadid = get_thread_num();

                for(j = 0; j < NTask; j++)
                  Thread[threadid].Exportflag[j] = -1;

                while(1)
                  {
                    if(Thread[threadid].Nexport >= (MaxNexport - (NTask - 1)) ||
              Thread[threadid].NexportNodes >= (MaxNexportNodes - NTopleaves))
                      break;

    #pragma omp atomic capture
                    idx = NextParticle++;

                    if(idx >= NActive)
                      break;

                    i = List[idx];

            tree_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
                  }
              }

              /* do all necessary bookkeeping and the data exchange */
              generic_exchange();


              /* now do the particles that were sent to us */
              int count = 0;
    #pragma omp parallel private(i)
              {
                int threadid = get_thread_num();

                while(1)
                  {
    #pragma omp atomic capture
                    i = count++;

                    if(i >= Nimport)
                      break;

                    tree_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
                  }
              }

              /* pick up results and process them */
              generic_finalize();

              /* check whether we are done */
              if(NextParticle >= NActive)
                ndone_flag = 1;
              else
                ndone_flag = 0;

              MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            }
          while(ndone < NTask);
    }


This assumes that the indices of the particles that you want to process are stored in a list called 
``List``, and that there are ``Nactive`` of them. The above is also enabled for use with OpenMP. 
(Note that OpenMP is not yet functional throughout the code, but it soon will be and then the above 
structures will deal with the complications involved in the parallel tree walks.) Almost all the 
statements should be kept as they are seen above, but the code to select the set of active particles 
can and should be changed as needed.


**[5]** Finally, we have to write the function that does the tree for a local or an imported particle. 
This is also where the non-trivial core of the computation happens. This function has the following 
generic template::

    static void tree_evaluate(int target, int mode, int threadid)
    {
      struct data_in local, *in;
      struct data_out out;

      if(mode == MODE_LOCAL_PARTICLES)
        {
          particle2in(&local, target, 0);
          in = &local;
          numnodes = 1;
          firstnode = NULL;
        }
      else
        {
          in = &DataGet[target];
          generic_get_numnodes(target, &numnodes, &firstnode);
        }


      /*  your code goes here. The input particle
       *  is stored in "in", the result should go to "out"
       */


      if(mode == MODE_LOCAL_PARTICLES)
        out2particle(&out, target, MODE_LOCAL_PARTICLES);
      else
        DataResult[target] = out;
    }


Again, this function should be kept as it is here. In the middle where the comment appears, your code 
should be added. It takes the information provided in the structure pointed to by ``in``, and fills 
in the result of the tree walks into the structure pointed to by ``out``. The starting nodes of the 
tree-branches that need to be walked are stored in the array ``firstnode[]``, and there ``numnodes`` 
of these. The tree walk can either be done by "hand" if needed, or by calling existing tree walk 
functions. For example, to get all gas cells within a radius ``Hsml``, you can call the function 
``ngb_treefind_variable_threads()``. The code for a SPH-like density estimate placed at the middle 
of the above function could then look something like::

    int nfound = ngb_treefind_variable_threads(in->Pos, in->Hsml, target, mode, threadid, numnodes, firstnode);

    out.Density = 0;

    for(n = 0; n < nfound; n++)
    {
        j = Thread[threadid].Ngblist[n];

        dx = in->Pos[0] - P[j].Pos[0];
        dy = in->Pos[1] - P[j].Pos[1];
        dz = in->Pos[2] - P[j].Pos[2];

        r2 = dx * dx + dy * dy + dz * dz;

        r = sqrt(r2);

        u = r / in->Hsml;

        if(u < 0.5)
          wk = hinv3 * (KERNEL_COEFF_1 + KERNEL_COEFF_2 * (u - 1) * u * u);
        else
          wk = hinv3 * KERNEL_COEFF_5 * (1.0 - u) * (1.0 - u) * (1.0 - u);

        out.Density += P[j].Mass * wk;
    }

   
**[6]** In some cases, there are no results returned from the tree walk. Then the part ::

    generic_exchange();

    /* ... secondary loop with MODE_IMPORTED_PARTICLES... */

    generic_finalize();

in the main routine can in principle be replaced by a single call of::

    generic_free_remaining_buffers();

Throughout the code, there are also a large number of examples to look at, but if anything is 
unclear, please ask.
