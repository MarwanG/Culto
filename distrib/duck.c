#include "duck.h"


t_int *duck_tilde_perform (t_int *w){
  return 0;
}

void duck_tilde_dsp(t_duck_tilde *x,t_signal **sp){
  dsp_add(duck_tilde_perform, 5, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_n); 
}


//DONE
void duck_tilde_free(t_duck_tilde *x){
  inlet_free(x->x_in2);
  outlet_free(x->x_out);
}


//MIGHT NEED SOMETHING FOR TSAMPLE
void *duck_tilde_new(void){
  t_duck_tilde *m;
  m =(t_duck_tilde *)pd_new(duck_tilde_class);
  m->x_in2 = inlet_new(&m->x_obj,&m->x_obj.ob_pd,&s_signal,&s_signal);
  m->x_out = outlet_new(&m->x_obj,&s_signal);
  return (void *)m;

}

//DONE
void duck_tilde_setup(void){
   duck_tilde_class = class_new(gensym("duck~"),
				(t_newmethod)duck_tilde_new, (t_method)duck_tilde_free
                              , sizeof(t_duck_tilde),
				 CLASS_DEFAULT, 0 );
   class_addmethod(duck_tilde_class, (t_method)duck_tilde_dsp,gensym("dsp"),0);
   CLASS_MAINSIGNALIN(duck_tilde_class, t_duck_tilde, f);  
}
