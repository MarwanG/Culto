#include "myftt~.h"


t_int* myftt_tilde_perform(t_int *w){
}


void myfft_tilde_dsp(t_myftt_tilde *x,t_signal **sp){


}

void myfft_tilde_free(t_myftt_tilde *x){
  outlet_free(x->x_out);
}


void *myfft_tilde_new(void){
  t_myftt_tilde *t;
  t = (t_myftt_tilde*)pd_new(myftt_tilde_class);
  t->x_out = outlet(&t->x_obj,&s_signal);
  return (void*)t;
}


void myfft_tilde_setup(void){
  myfft_tilde_class = class_new(gensym("myftt~",(t_newmethod)myfft_tilde_new,
				       (t_method)myfft_tilde_free,
				       sizeof(t_myfft_tilde),
				       CLASS_DEFAULT,0));
  class_addmethod(myfft_tilde_class,(t_method)duck_tilde_dsp,gensym("dsp"),0);

  CLASS_MAINSIGNALIN(myfft_tilde_class,t_duck_tilde,f);
 
}
