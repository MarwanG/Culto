#include <math.h>
#include <stdlib.h>
#include "myfft~.h"
#define BUFFER_LEN 4096

t_int* myftt_tilde_perform(t_int *w){
  int i;
  t_myfft_tilde *x = (t_myfft_tilde *) w[1];
  
  t_sample* vec1 = (t_sample*)w[2];
  int * buffer = x->buffer;
  t_sample* output = (t_sample*)w[4];
  for(i=0;i<w[5];i++){
    buffer[x->cpt+i] = vec1[i];
  }
  x->cpt = x->cpt + w[5];
  if(x->cpt >= BUFFER_LEN){
    /*Applic un fenetre de Hamming*/
    int T = BUFFER_LEN; /*Je suis pas sur c pas marque et on ligne c le length*/
    for(i=0;i<BUFFER_LEN;i++){
      if(buffer[i] > BUFFER_LEN || buffer[i] < 0){
	    buffer[i] = 0.54 - (0.46 * cos(2 * PI *(buffer[i]/T)));
      }else{
	    buffer[i] = 0;
      }
    }
    for(i=0;i<BUFFER_LEN;i++){
      output[i] = buffer[i];
    }
  }
  return w+6;
}


void myfft_tilde_dsp(t_myfft_tilde *x,t_signal **sp){
 dsp_add(myftt_tilde_perform, 5, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_n); 
}

void myfft_tilde_free(t_myfft_tilde *x){
  outlet_free(x->x_out);
}


void *myfft_tilde_new(void){
  t_myfft_tilde *t;
  t = (t_myfft_tilde*)pd_new(myfft_tilde_class);
  t->x_out = outlet_new(&t->x_obj,&s_signal);
  printf("kaka");
  t->buffer = (int*)malloc(BUFFER_LEN * sizeof (int));
  t->cpt = 0;
  return (void*)t;
}


void myfft_tilde_setup(void){
  myfft_tilde_class = class_new(gensym("myfft~"),(t_newmethod)myfft_tilde_new,
				       (t_method)myfft_tilde_free,
				       sizeof(t_myfft_tilde),
				       CLASS_DEFAULT,0);
  class_addmethod(myfft_tilde_class,(t_method)myfft_tilde_dsp,gensym("dsp"),0);

  CLASS_MAINSIGNALIN(myfft_tilde_class,t_myfft_tilde,f);
 
}


    
