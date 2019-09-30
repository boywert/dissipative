
GFM_CHEMTAGS
============

A tagging/tracking scheme which can follow how much metal mass is separately ejected in SNIa, SNII, and AGB processes, 
as well as total mass ejected from NS-NS mergers, and iron (Fe) mass from SNIa and SNII separately. In all cases, 
evolving stellar populations inject mass tracers into their surrounding gas environment.
This information is output in a new variable ``GFM_Chemtags`` for both gas and star particles. Note: for stars, 
all values are simply inherited at the time of formation from the gas cell from which the star was born. They do not 
then evolve or change in any way (i.e. no self-enrichment), so these values describe the 'inherited' wind/SN/NSNS material 
from the gas. 


======  ============================  ===================================================================================
Entry   Units                         Description
======  ============================  ===================================================================================
SNIa    dimensionless mass ratio      The total metal mass ejected by Type Ia supernovae.
SNII    dimensionless mass ratio      The total metal mass ejected by Type II supernovae.
AGB     dimensionless mass ratio      The total metal mass ejected by stellar winds, which is dominated by AGB stars.
NSNS    same as `NSNS_MassPerEvent`   The total mass ejected from NS-NS merger events (see below).
FeSNIA  dimensionless mass ratio      The total iron ejected by Type Ia supernovae alone.
FeSNII  dimensionless mass ratio      The total iron ejected by Type II supernovae alone.
======  ============================  ===================================================================================

Depending on Config options, ``GFM_ChemTags`` may contain 3, 4, 5 or 6 entries per gas/star. In these cases, the entries 
(and their order) are::

  3: SNIa, SNII, AGB
  4: SNIa, SNII, AGB, NSNS
  5: SNIa, SNII, AGB, FeSNIA, FeSNII
  6: SNIa, SNII, AGB, NSNS, FeSNIA, FeSNII (as in TNG)

The units of all these fields except for NSNS are the same as GFM_Metals: dimensionless mass ratios. 
If you sum up all the elements of ``GFM_Metals`` heavier than Helium, you recover the sum of the three tags SNIa+SNII+AGB. 
Likewise, the Fe entry of ``GFM_Metals`` roughly equals the sum of FeSNIa+FeSNII, modulo the small amount of iron consumed 
(i.e. negative contribution) by AGB winds.


Usage
-----

To use, enable ``GFM_CHEMTAGS`` in the Config file, and add ``256`` to ``GFM_OUTPUT_MASK`` to save the tags in the snapshot files.


Additional Parameters
---------------------

* ``NSNS_MassTransferOn`` Master switch. Always set to 1.
* ``NSNS_MassPerEvent`` The mass ejected per NS-NS merger event (e.g. 0.05 solar masses). The units of the corresponding 
  output entry are the same as whatever units are chosen here for this parameter. However, this parameter could be e.g. solar 
  masses times some large factor :math:`\alpha`, in which case this output field should be normalized by :math:`\alpha` to 
  recover physical solar masses. Likewise, for any preferred physical value of the mass ejected per NS-NS merger, this can 
  be enforced in post-processing by multiplying this output field by (MyPreferred_NSNS_MassPerEvent / NSNS_MassPerEvent). 
  In practice, this boost factor approach was (or should be) used, as it was found that such a small value led to numerical 
  truncation issues during advection.
  *Recommended value: 5000*. (0.05 solar masses from Shen+ 2015 times a boost factor of 100000)
* ``NSNS_Rate_TAU`` The timescale :math:`\tau` for the powerlaw DTD for NS-NS merger events.
  *Recommended/TNG value: 0.1*. (from Shen+ 2015) (0.04 would be the same as SNIa)
* ``NSNS_per_SNIa`` The number of NS-NS mergers per SNIa.
  *Recommended/TNG value: 1e-3*. (see Shen+ 2015, references therein, and JN/ERR)


Additional Config.sh Options
----------------------------

GFM_CHEMTAGS
  Master switch, to track extra species in any manner (with/without split Fe, with/without NSNS tracking, etc).

GFM_SPLITFE
  This toggle will split the Fe abudance, and track iron produced by SNIa and SNII separately.

GFM_RPROCESS
  This toggle enables tracking total mass ejecta from NS-NS mergers. NS-NS merger events are modeled 
  stochastically (i.e. no fractional events) with a DTD scheme similar to that used for SNIa, except with 
  a different :math:`\tau` value. If enabled, then must also specify the parameters ``NSNS_MassPerEvent``, 
  ``NSNS_Rate_TAU``, and ``NSNS_per_SNIa``. Note: the DTD powerlaw index is hardcoded like the SNIa power law index.

GFM_SPLITFE_ADDINAGB
  Status unknown, do not use. "add in the AGB iron half-half on the two iron SNIa/SNII tags such that the sum of 
  them should be equal to the total iron".

Authors
-------

  * Jill Naiman (jill.naiman@cfa.harvard.edu)
  * Volker Springel (volker.springel@h-its.org)


Usage Policy and Citation
-------------------------

Please contact the authors before using this code for a new project. 

Papers to cite:

  * Naiman et al. (in prep) (TNG first use of scheme)
  * Pillepich et al. (in prep) (TNG methods2)
