#include "sms~.h"
#include <stdlib.h>
#define BUFFER_LEN 4096

t_int *sms_tilde_perform (t_int *w){
    t_sms_tilde* x = (t_sms_tilde*) w[1];
    t_inlet* in_srcHarm = (t_inlet*) w[2];
    t_inlet* in_srcMod = (t_inlet*) w[3];
    t_inlet* in_synth = (t_inlet*) w[4];
    t_inlet* in_msg = (t_inlet*) w[5];
    t_outlet* out = (t_outlet*) w[6];
    int vecSize = (int) w[7];

    if (in_msg[1]->info == 1){ // if bypass
        
    }

    

    printf("coucou\n");
    return 0;
}

void sms_tilde_dsp(t_sms_tilde *x,t_signal **sp){
    dsp_add(sms_tilde_perform, 7, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, sp[4]->s_vec, sp[0]->s_n); 

}

void sms_tilde_free(t_sms_tilde *x){
	inlet_free(x->in1);
	inlet_free(x->in2);
	inlet_free(x->in3);
	inlet_free(x->in4);
	outlet_free(x->out);
	free(x->buffer);
}

void *sms_tilde_new(t_symbol *s, int argc, t_atom *argv){
	t_sms_tilde *m;
	m = (t_sms_tilde*)pd_new(sms_tilde_class);
	m->in1 = inlet_new(&m->x_obj,&m->x_obj.ob_pd,&s_signal,&s_signal);
	m->in2 = inlet_new(&m->x_obj,&m->x_obj.ob_pd,&s_signal,&s_signal);
 	m->in3 = floatinlet_new(&m->x_obj,&m->seuil);
 	m->in4 = floatinlet_new(&m->x_obj,&m->info);
 	m->out = outlet_new(&m->x_obj,&s_signal);

 	if(argc < 2){
 		m->buffer = (t_sample*)malloc(BUFFER_LEN * sizeof (t_sample));
 	}else{
 		m->buffer = (t_sample*)malloc(atom_getint(argv) * sizeof (t_sample));
 	}

 	return (void*) m;
}

void sms_tilde_setup(void){
	sms_tilde_class = class_new(gensym("sms~"),(t_newmethod)sms_tilde_new
		,(t_method)sms_tilde_free,sizeof(t_sms_tilde),
		CLASS_DEFAULT,0);
	class_addmethod(sms_tilde_class, (t_method)sms_tilde_dsp,gensym("dsp"),0);
   	CLASS_MAINSIGNALIN(sms_tilde_class, t_sms_tilde, f);  

}
