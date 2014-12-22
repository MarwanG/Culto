#include "myftt~.h"
#include <math.h>


t_int* myftt_tilde_perform(t_int *w){
  int i;
  t_sample* vec1 = (t_sample*)w[2];
  int * buffer = w[0]->buffer;
  for(i=0;i<w[5];i++){
    buffer[w[0]->cpt+i] = vec1[i];
  }
  cpt = cpt + w[5];
  if(cpt >= 2048){
    /*Applic un fenetre de Hamming*/
    int T = 2048; /*Je suis pas sur c pas marque et on ligne c le length*/
    for(i=0;i<2048;i++){
      if(buffer[i] > 2048 || buffer[i] < 0){
	buffer[i] = 0.54 - (0.46 * cos(2 * PI *(buffer[i]/T)));
      }else{
	buffer[i] = 0;
      }
    }
  }
  return w+6;
}


void myfft_tilde_dsp(t_myftt_tilde *x,t_signal **sp){
 dsp_add(myftt_tilde_perform, 5, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_n); 
}

void myfft_tilde_free(t_myftt_tilde *x){
  outlet_free(x->x_out);
}


void *myfft_tilde_new(void){
  t_myftt_tilde *t;
  t = (t_myftt_tilde*)pd_new(myftt_tilde_class);
  t->x_out = outlet(&t->x_obj,&s_signal);
  t->buffer = malloc(2048 * sizeof int);
  t->cpt = 0;
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


    
