/*!
 * \copyright   This file is part of the AREPO code developed by Volker Springel.
 * \copyright   Copyright (C) 2013  by Volker Springel (volker.springel@h-its.org)
 * \copyright   and contributing authors.
 *  
 * \file        src/fm_star_formation/early_stellar_feedback_proto.h
 * \date        MM/YYYY
 * \author     
 * \brief        
 * \details     
 * 
 * 
 * \par Major modifications and contributions:
 * 
 * - DD.MM.YYYY Description
 */

#ifndef FM_EARLY_STELLAR_FEED_PROTO_H
#define FM_EARLY_STELLAR_FEED_PROTO_H

int is_doing_early_feedback(int i);

MyDouble compute_EarlyFeedback_energy(MyDouble age_star_in_gyr, MyFloat mass, MyFloat dt);

void do_early_stellar_feedback(void);
int early_stellar_feedback_evaluate(int target, int mode, int thread_id);
#endif
