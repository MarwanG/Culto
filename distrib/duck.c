#include "duck.h"


t_int *duck_tilde_perform (t_int *w){

}

void duck_tilde_dsp(t_duck_tilde *x,t_signal **sp){

}

void duck_tilde_free(t_pan_tilde *x){
  inlet_free(x->x_in2);
  outlet_free(x->x_out);
}

void *duck_tilde_new(void){
  t_duck_tilde *m;
  m =(t_duck_tilde *)pd_new(duck_tilde_class);
  m->x_in2 = inlet_new(m,&m->x_obj.ob_pd,&s_signal,&s_signal);
  m->x_out = outlet_new(m,&s_signal);
  /* t_sample WTF to do*/
  return (void *)m;

}

void duck_tilde_setup(void){
   duck_tilde_class = class_new(gensym("duck_tilde"),
				(t_newmethod)duck_tilde_new, duck_tilde_free
                              , sizeof(t_duck_tilde),
				 CLASS_DEFAULT, 0 );
   class_addmethod(duck_tilde_class, duck_tilde_dsp,gensym("dsp"),0);
   CLASS_MAINSIGNALIN(duck_tilde_class, t_duck_tilde, f);  
}
