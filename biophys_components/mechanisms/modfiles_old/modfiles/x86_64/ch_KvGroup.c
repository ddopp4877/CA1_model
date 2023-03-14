/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__ch_KvGroup
#define _nrn_initial _nrn_initial__ch_KvGroup
#define nrn_cur _nrn_cur__ch_KvGroup
#define _nrn_current _nrn_current__ch_KvGroup
#define nrn_jacob _nrn_jacob__ch_KvGroup
#define nrn_state _nrn_state__ch_KvGroup
#define _net_receive _net_receive__ch_KvGroup 
#define _f_trates _f_trates__ch_KvGroup 
#define rates rates__ch_KvGroup 
#define states states__ch_KvGroup 
#define trates trates__ch_KvGroup 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gmax _p[0]
#define g _p[1]
#define ik _p[2]
#define ninf _p[3]
#define ntau _p[4]
#define myi _p[5]
#define n _p[6]
#define ek _p[7]
#define Dn _p[8]
#define nexp _p[9]
#define _g _p[10]
#define _ion_ek	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_rates(void);
 static void _hoc_states(void);
 static void _hoc_trates(void);
 static void _hoc_vtrap(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_ch_KvGroup", _hoc_setdata,
 "rates_ch_KvGroup", _hoc_rates,
 "states_ch_KvGroup", _hoc_states,
 "trates_ch_KvGroup", _hoc_trates,
 "vtrap_ch_KvGroup", _hoc_vtrap,
 0, 0
};
#define vtrap vtrap_ch_KvGroup
 extern double vtrap( double , double );
 /* declare global and static user variables */
#define usetable usetable_ch_KvGroup
 double usetable = 1;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "usetable_ch_KvGroup", 0, 1,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "gmax_ch_KvGroup", "mho/cm2",
 "g_ch_KvGroup", "mho/cm2",
 "ik_ch_KvGroup", "mA/cm2",
 "ntau_ch_KvGroup", "ms",
 "myi_ch_KvGroup", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double n0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "usetable_ch_KvGroup", &usetable_ch_KvGroup,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"ch_KvGroup",
 "gmax_ch_KvGroup",
 0,
 "g_ch_KvGroup",
 "ik_ch_KvGroup",
 "ninf_ch_KvGroup",
 "ntau_ch_KvGroup",
 "myi_ch_KvGroup",
 0,
 "n_ch_KvGroup",
 0,
 0};
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 11, _prop);
 	/*initialize range parameters*/
 	gmax = 0;
 	_prop->param = _p;
 	_prop->param_size = 11;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _ch_KvGroup_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("k", 1.0);
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 11, 3);
  hoc_register_dparam_semantics(_mechtype, 0, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "k_ion");
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 ch_KvGroup /home/dpd4k4/CA1_model/biophys_components/mechanisms/modfiles/ch_KvGroup.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
 static double FARADAY = 96520.0;
 static double R = 8.3134;
 static double _zq10 ;
 static double *_t_ninf;
 static double *_t_nexp;
 static double *_t_ntau;
static int _reset;
static char *modelname = "Slower potassium channels (voltage dependent)";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int _f_trates(double);
static int rates(double);
static int states();
static int trates(double);
 static void _n_trates(double);
 
/*VERBATIM*/
#include <stdlib.h> /* 	Include this library so that the following
						(innocuous) warning does not appear:
						 In function '_thread_cleanup':
						 warning: incompatible implicit declaration of 
						          built-in function 'free'  */
 
static int  states (  ) {
   trates ( _threadargscomma_ v ) ;
   n = n + nexp * ( ninf - n ) ;
    return 0; }
 
static void _hoc_states(void) {
  double _r;
   _r = 1.;
 states (  );
 hoc_retpushx(_r);
}
 
static int  rates (  double _lv ) {
   double _lalpha , _lbeta , _lsum ;
 _zq10 = pow( 3.0 , ( ( celsius - 34.0 ) / 10.0 ) ) ;
   _lalpha = 0.0189324 * vtrap ( _threadargscomma_ - ( _lv - 4.18371 ) , 6.42606 ) ;
   _lbeta = 0.015857 * exp ( - _lv / 25.4834 ) ;
   _lsum = _lalpha + _lbeta ;
   ntau = 1.0 / _lsum ;
   ninf = _lalpha / _lsum ;
    return 0; }
 
static void _hoc_rates(void) {
  double _r;
   _r = 1.;
 rates (  *getarg(1) );
 hoc_retpushx(_r);
}
 static double _mfac_trates, _tmin_trates;
 static void _check_trates();
 static void _check_trates() {
  static int _maktable=1; int _i, _j, _ix = 0;
  double _xi, _tmax;
  static double _sav_dt;
  static double _sav_celsius;
  if (!usetable) {return;}
  if (_sav_dt != dt) { _maktable = 1;}
  if (_sav_celsius != celsius) { _maktable = 1;}
  if (_maktable) { double _x, _dx; _maktable=0;
   _tmin_trates =  - 100.0 ;
   _tmax =  100.0 ;
   _dx = (_tmax - _tmin_trates)/200.; _mfac_trates = 1./_dx;
   for (_i=0, _x=_tmin_trates; _i < 201; _x += _dx, _i++) {
    _f_trates(_x);
    _t_ninf[_i] = ninf;
    _t_nexp[_i] = nexp;
    _t_ntau[_i] = ntau;
   }
   _sav_dt = dt;
   _sav_celsius = celsius;
  }
 }

 static int trates(double _lv){ _check_trates();
 _n_trates(_lv);
 return 0;
 }

 static void _n_trates(double _lv){ int _i, _j;
 double _xi, _theta;
 if (!usetable) {
 _f_trates(_lv); return; 
}
 _xi = _mfac_trates * (_lv - _tmin_trates);
 if (isnan(_xi)) {
  ninf = _xi;
  nexp = _xi;
  ntau = _xi;
  return;
 }
 if (_xi <= 0.) {
 ninf = _t_ninf[0];
 nexp = _t_nexp[0];
 ntau = _t_ntau[0];
 return; }
 if (_xi >= 200.) {
 ninf = _t_ninf[200];
 nexp = _t_nexp[200];
 ntau = _t_ntau[200];
 return; }
 _i = (int) _xi;
 _theta = _xi - (double)_i;
 ninf = _t_ninf[_i] + _theta*(_t_ninf[_i+1] - _t_ninf[_i]);
 nexp = _t_nexp[_i] + _theta*(_t_nexp[_i+1] - _t_nexp[_i]);
 ntau = _t_ntau[_i] + _theta*(_t_ntau[_i+1] - _t_ntau[_i]);
 }

 
static int  _f_trates (  double _lv ) {
   double _ltinc ;
 rates ( _threadargscomma_ _lv ) ;
   _ltinc = - dt * _zq10 ;
   nexp = 1.0 - exp ( _ltinc / ntau ) ;
    return 0; }
 
static void _hoc_trates(void) {
  double _r;
    _r = 1.;
 trates (  *getarg(1) );
 hoc_retpushx(_r);
}
 
double vtrap (  double _lx , double _ly ) {
   double _lvtrap;
 if ( fabs ( _lx / _ly ) < 1e-6 ) {
     _lvtrap = _ly * ( 1.0 - _lx / _ly / 2.0 ) ;
     }
   else {
     _lvtrap = _lx / ( exp ( _lx / _ly ) - 1.0 ) ;
     }
   
return _lvtrap;
 }
 
static void _hoc_vtrap(void) {
  double _r;
   _r =  vtrap (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ hoc_execerror("ch_KvGroup", "cannot be used with CVODE"); return 0;}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_k_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  n = n0;
 {
   trates ( _threadargscomma_ v ) ;
   n = ninf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ek = _ion_ek;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   g = gmax * n ;
   ik = g * ( v - ek ) ;
   myi = ik ;
   }
 _current += ik;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ek = _ion_ek;
 _g = _nrn_current(_v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ek = _ion_ek;
 { error =  states();
 if(error){fprintf(stderr,"at line 70 in file ch_KvGroup.mod:\n	SOLVE states\n"); nrn_complain(_p); abort_run(error);}
 } }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
   _t_ninf = makevector(201*sizeof(double));
   _t_nexp = makevector(201*sizeof(double));
   _t_ntau = makevector(201*sizeof(double));
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/home/dpd4k4/CA1_model/biophys_components/mechanisms/modfiles/ch_KvGroup.mod";
static const char* nmodl_file_text = 
  "TITLE Slower potassium channels (voltage dependent)\n"
  "\n"
  "COMMENT\n"
  "Slower potassium channels (voltage dependent)\n"
  "\n"
  "Ions: k\n"
  "\n"
  "Style: quasi-ohmic\n"
  "\n"
  "From: Modified from Yuen and Durand, 1991 (squid axon) by Aradi and Holmes (2002)\n"
  "\n"
  "Updates:\n"
  "2014 December (Marianne Bezaire): documented\n"
  "2014 February (Marianne Bezaire): modified to fit the experimental data from Lu and Jonas\n"
  "? ? ?: shifted by 65 mV as compared to the Aradi & Holmes paper, unsure why\n"
  "ENDCOMMENT\n"
  "\n"
  "\n"
  "VERBATIM\n"
  "#include <stdlib.h> /* 	Include this library so that the following\n"
  "						(innocuous) warning does not appear:\n"
  "						 In function '_thread_cleanup':\n"
  "						 warning: incompatible implicit declaration of \n"
  "						          built-in function 'free'  */\n"
  "ENDVERBATIM\n"
  " \n"
  "UNITS {\n"
  "	(mA) =(milliamp)\n"
  "	(mV) =(millivolt)\n"
  "	(uF) = (microfarad)\n"
  "	(molar) = (1/liter)\n"
  "	(nA) = (nanoamp)\n"
  "	(mM) = (millimolar)\n"
  "	(um) = (micron)\n"
  "	FARADAY = 96520 (coul)\n"
  "	R = 8.3134	(joule/degC)\n"
  "}\n"
  " \n"
  "NEURON { \n"
  "	SUFFIX ch_KvGroup\n"
  "	USEION k READ ek WRITE ik VALENCE 1\n"
  "	RANGE g, gmax, ninf, ntau, ik\n"
  "	RANGE myi\n"
  "	THREADSAFE\n"
  "}\n"
  " \n"
  "PARAMETER {\n"
  "	v (mV) \n"
  "	celsius (degC) : temperature - set in hoc; default is 6.3\n"
  "	dt (ms) \n"
  "\n"
  "	ek  (mV)\n"
  "	gmax (mho/cm2)\n"
  "}\n"
  " \n"
  "STATE {\n"
  "	n	\n"
  "}\n"
  " \n"
  "ASSIGNED {		     \n"
  "	g (mho/cm2)\n"
  "	ik (mA/cm2)\n"
  "	ninf\n"
  "	ntau (ms)\n"
  "	nexp\n"
  "	myi (mA/cm2)\n"
  "} \n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states\n"
  "	g = gmax*n\n"
  "	ik = g*(v-ek)\n"
  "	myi =  ik\n"
  "}\n"
  " \n"
  "UNITSOFF\n"
  " \n"
  "INITIAL {\n"
  "	trates(v)\n"
  "\n"
  "	n = ninf\n"
  "}\n"
  "\n"
  "PROCEDURE states() {	:Computes state variables m, h, and n \n"
  "	trates(v)	:      at the current v and dt.       \n"
  "	n = n + nexp*(ninf-n)\n"
  "}\n"
  " \n"
  "LOCAL q10\n"
  "PROCEDURE rates(v) {  :Computes rate and other constants at current v.\n"
  "                      :Call once from HOC to initialize inf at resting v.\n"
  "	LOCAL  alpha, beta, sum\n"
  "	q10 = 3^((celsius - 34)/10)\n"
  "		\n"
  "	:\"nf\" fKDR activation system\n"
  "	:alpha = a*(-(v+b))/(exp(-(v+b)/c)-1) : Lien and Jonas, 2003\n"
  "	alpha = 0.0189324*vtrap(-(v-4.18371),6.42606) : Lien and Jonas, 2003\n"
  "	beta = 0.015857*exp(-v/25.4834) : Lien and Jonas, 2003\n"
  "\n"
  "	sum = alpha+beta        \n"
  "	ntau = 1/sum\n"
  "	ninf = alpha/sum	\n"
  "}\n"
  " \n"
  "PROCEDURE trates(v) {  :Computes rate and other constants at current v.\n"
  "                      :Call once from HOC to initialize inf at resting v.\n"
  "	LOCAL tinc\n"
  "	TABLE ninf, nexp, ntau\n"
  "	DEPEND dt, celsius\n"
  "	FROM -100 TO 100 WITH 200\n"
  "						   \n"
  "	rates(v)	: not consistently executed from here if usetable_hh == 1\n"
  "	: so don't expect the tau values to be tracking along with\n"
  "	: the inf values in hoc\n"
  "\n"
  "	tinc = -dt * q10\n"
  "	nexp = 1 - exp(tinc/ntau)\n"
  "}\n"
  " \n"
  "FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.\n"
  "        if (fabs(x/y) < 1e-6) {\n"
  "                vtrap = y*(1 - x/y/2)\n"
  "        }else{  \n"
  "                vtrap = x/(exp(x/y) - 1)\n"
  "        }\n"
  "}\n"
  " \n"
  "UNITSON\n"
  "\n"
  ;
#endif