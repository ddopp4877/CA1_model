#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _bg2pyr_reg(void);
extern void _ch_CavL_reg(void);
extern void _ch_CavN_reg(void);
extern void _ch_HCN_reg(void);
extern void _ch_HCNolm_reg(void);
extern void _ch_HCNp_reg(void);
extern void _ch_KCaS_reg(void);
extern void _ch_Kdrfast_reg(void);
extern void _ch_Kdrfastngf_reg(void);
extern void _ch_Kdrp_reg(void);
extern void _ch_Kdrslow_reg(void);
extern void _ch_KvAdistp_reg(void);
extern void _ch_KvA_reg(void);
extern void _ch_KvAngf_reg(void);
extern void _ch_KvAolm_reg(void);
extern void _ch_KvAproxp_reg(void);
extern void _ch_KvCaB_reg(void);
extern void _ch_KvGroup_reg(void);
extern void _ch_KvM_reg(void);
extern void _ch_leak_reg(void);
extern void _chn2pyr_reg(void);
extern void _ch_Navaxonp_reg(void);
extern void _ch_Navbis_reg(void);
extern void _ch_Navcck_reg(void);
extern void _ch_Nav_reg(void);
extern void _ch_Navngf_reg(void);
extern void _ch_Navp_reg(void);
extern void _ExpGABAab_reg(void);
extern void _iconc_Ca_reg(void);
extern void _int2int_reg(void);
extern void _int2pyr_reg(void);
extern void _kv_reg(void);
extern void _MyExp2Sid_reg(void);
extern void _MyExp2Sidnw_reg(void);
extern void _mynetstim_reg(void);
extern void _na12_reg(void);
extern void _na16_reg(void);
extern void _positionfcns_reg(void);
extern void _pyr2int_reg(void);
extern void _pyr2pyr_reg(void);
extern void _SIN_reg(void);
extern void _vecevent_reg(void);
extern void _xtra_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," \"bg2pyr.mod\"");
    fprintf(stderr," \"ch_CavL.mod\"");
    fprintf(stderr," \"ch_CavN.mod\"");
    fprintf(stderr," \"ch_HCN.mod\"");
    fprintf(stderr," \"ch_HCNolm.mod\"");
    fprintf(stderr," \"ch_HCNp.mod\"");
    fprintf(stderr," \"ch_KCaS.mod\"");
    fprintf(stderr," \"ch_Kdrfast.mod\"");
    fprintf(stderr," \"ch_Kdrfastngf.mod\"");
    fprintf(stderr," \"ch_Kdrp.mod\"");
    fprintf(stderr," \"ch_Kdrslow.mod\"");
    fprintf(stderr," \"ch_KvAdistp.mod\"");
    fprintf(stderr," \"ch_KvA.mod\"");
    fprintf(stderr," \"ch_KvAngf.mod\"");
    fprintf(stderr," \"ch_KvAolm.mod\"");
    fprintf(stderr," \"ch_KvAproxp.mod\"");
    fprintf(stderr," \"ch_KvCaB.mod\"");
    fprintf(stderr," \"ch_KvGroup.mod\"");
    fprintf(stderr," \"ch_KvM.mod\"");
    fprintf(stderr," \"ch_leak.mod\"");
    fprintf(stderr," \"chn2pyr.mod\"");
    fprintf(stderr," \"ch_Navaxonp.mod\"");
    fprintf(stderr," \"ch_Navbis.mod\"");
    fprintf(stderr," \"ch_Navcck.mod\"");
    fprintf(stderr," \"ch_Nav.mod\"");
    fprintf(stderr," \"ch_Navngf.mod\"");
    fprintf(stderr," \"ch_Navp.mod\"");
    fprintf(stderr," \"ExpGABAab.mod\"");
    fprintf(stderr," \"iconc_Ca.mod\"");
    fprintf(stderr," \"int2int.mod\"");
    fprintf(stderr," \"int2pyr.mod\"");
    fprintf(stderr," \"kv.mod\"");
    fprintf(stderr," \"MyExp2Sid.mod\"");
    fprintf(stderr," \"MyExp2Sidnw.mod\"");
    fprintf(stderr," \"mynetstim.mod\"");
    fprintf(stderr," \"na12.mod\"");
    fprintf(stderr," \"na16.mod\"");
    fprintf(stderr," \"positionfcns.mod\"");
    fprintf(stderr," \"pyr2int.mod\"");
    fprintf(stderr," \"pyr2pyr.mod\"");
    fprintf(stderr," \"SIN.mod\"");
    fprintf(stderr," \"vecevent.mod\"");
    fprintf(stderr," \"xtra.mod\"");
    fprintf(stderr, "\n");
  }
  _bg2pyr_reg();
  _ch_CavL_reg();
  _ch_CavN_reg();
  _ch_HCN_reg();
  _ch_HCNolm_reg();
  _ch_HCNp_reg();
  _ch_KCaS_reg();
  _ch_Kdrfast_reg();
  _ch_Kdrfastngf_reg();
  _ch_Kdrp_reg();
  _ch_Kdrslow_reg();
  _ch_KvAdistp_reg();
  _ch_KvA_reg();
  _ch_KvAngf_reg();
  _ch_KvAolm_reg();
  _ch_KvAproxp_reg();
  _ch_KvCaB_reg();
  _ch_KvGroup_reg();
  _ch_KvM_reg();
  _ch_leak_reg();
  _chn2pyr_reg();
  _ch_Navaxonp_reg();
  _ch_Navbis_reg();
  _ch_Navcck_reg();
  _ch_Nav_reg();
  _ch_Navngf_reg();
  _ch_Navp_reg();
  _ExpGABAab_reg();
  _iconc_Ca_reg();
  _int2int_reg();
  _int2pyr_reg();
  _kv_reg();
  _MyExp2Sid_reg();
  _MyExp2Sidnw_reg();
  _mynetstim_reg();
  _na12_reg();
  _na16_reg();
  _positionfcns_reg();
  _pyr2int_reg();
  _pyr2pyr_reg();
  _SIN_reg();
  _vecevent_reg();
  _xtra_reg();
}

#if defined(__cplusplus)
}
#endif
