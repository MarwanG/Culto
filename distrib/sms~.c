#include "sms~.h"
#define BUFFER_LEN 4096

_int *sms_tilde_perform (t_int *w){

}

void sms_tilde_dsp(t_duck_tilde *x,t_signal **sp){

}

void sms_tilde_free(t_duck_tilde *x){
	inlet_free(x->in1);
	inlet_free(x->in2);
	inlet_free(x->in3);
	inlet_free(x->in4);
	outlet_free(x->x_out);
	free(x->buffer);
}

void *sms_tilde_new(t_symbol *s, int argc, t_atom *argv){
	t_sms_tilde *m;
	m = (t_sms_tilde*)pd_new(sms_tilde_class);
	m->x_in1 = inlet_new(&m->x_obj,&m->x_obj.ob_pd,&s_signal,&s_signal);
	m->x_in2 = inlet_new(&m->x_obj,&m->x_obj.ob_pd,&s_signal,&s_signal);
 	m->x_in3 = floatinlet_new(&m->x_obj,&m->seuil);
  	m->x_in4 = pointerinlet_new(&m->x_obj,&m->info);
 	m->x_out = outlet_new(&m->x_obj,&s_signal);

 	if(argc == 1){
 		m->buffer = (t_sample*)malloc(atom_getint(argv) * sizeof (t_sample));
 	}else{
 		m->buffer = (t_sample*)malloc(BUFFER_LEN * sizeof (t_sample));
 	}

 	return (void*m);
}

void sms_tilde_setup(void){
	sms_tilde_class = class_new(gensym("sms~"),(t_newmethod)sms_tilde_new
		,(t_method)sms_tilde_free,sizeof(t_sms_tilde),
		CLASS_DEFAULT,0);
	class_addmethod(sms_tilde_class, (t_method)sms_tilde_dsp,gensym("dsp"),0);
   	CLASS_MAINSIGNALIN(sms_tilde_class, t_sms_tilde, f);  

}