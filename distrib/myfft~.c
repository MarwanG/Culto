#include <math.h>
#include <stdlib.h>
#include "myfft~.h"
#define BUFFER_LEN 4096

t_int* myftt_tilde_perform(t_int *w){
  int i;
  t_myfft_tilde *x = (t_myfft_tilde *) (w[1]);
  int * buffer = (int*)x->buffer;
  printf("coucou1\n");
  printf("%d\n",x->cpt);

  
  t_sample* vec1 = (t_sample*) w[2];
  /*vec1 = (t_sample*) realloc(vec1, (2*BUFFER_LEN)*sizeof(t_sample));
  if (vec1 == NULL){
    perror("malloc");
    exit(EXIT_FAILURE);
  }*/
  t_sample* output = (t_sample*) w[3];
  /*output = (t_sample*) realloc(output, (2*BUFFER_LEN)*sizeof(t_sample));
  if (output == NULL){
    perror("malloc");
    exit(EXIT_FAILURE);
  }*/
  printf("coucou2\n");
  for(i=0;i<w[4];i++){
    buffer[x->cpt+i] = vec1[i];
  }
  if(x->cpt+w[4] >= BUFFER_LEN){
    /*Applic un fenetre de Hamming*/
    for(i=x->cpt;i<BUFFER_LEN;i++){
      if(buffer[i] < BUFFER_LEN || buffer[i] < 0){
	    buffer[i] = 0.54 - (0.46 * cos(2 * PI *(buffer[i]/BUFFER_LEN)));
      }else{
	    buffer[i] = 0;
      }
    }
    x->cpt = x->cpt + w[4];
    for(i=0;i<w[4];i++){
      output[i] = buffer[i];
    }
    int j=0;
    for(i=w[4];i<BUFFER_LEN;i++){
        buffer[j] = buffer[i];
        j++;
    }
    x->cpt = j;
    init_rdft(w[4], x->bitshuffle, x->weighting);
    rdft(w[4], 0, output, x->bitshuffle, x->weighting);
  }

  return w+5;
}


void myfft_tilde_dsp(t_myfft_tilde *x,t_signal **sp){
  printf("tilde_dsp\n");
  printf("point %p\n",(void*)x);
  dsp_add(myftt_tilde_perform, 4, x, sp[0]->s_vec, sp[1]->s_vec, sp[0]->s_n); 
  printf("tilde_dsp done\n");
}

void myfft_tilde_free(t_myfft_tilde *x){
  printf("tilde_free\n");
  outlet_free(x->x_out);
  free(x->buffer);
}


void *myfft_tilde_new(void){
  printf("tile_new\n");
  t_myfft_tilde *t;
  t = (t_myfft_tilde*)pd_new(myfft_tilde_class);
  t->x_out = outlet_new(&t->x_obj,&s_signal);
  printf("coucou new 1\n");
  t->buffer = (t_sample*)malloc(BUFFER_LEN * sizeof (t_sample));
  t->cpt = 0;
  printf("coucou new 2\n");
  return (void*)t;
}


void myfft_tilde_setup(void){
  myfft_tilde_class = class_new(gensym("myfft~"),(t_newmethod)myfft_tilde_new,
				       (t_method)myfft_tilde_free,
				       sizeof(t_myfft_tilde),
				       CLASS_DEFAULT,0);
  printf("tile_newd done\n");
  class_addmethod(myfft_tilde_class,(t_method)myfft_tilde_dsp,gensym("dsp"),0);
  printf("tile_dsp done\n");

  CLASS_MAINSIGNALIN(myfft_tilde_class,t_myfft_tilde,f);
 
  printf("tile_setup done\n");
}


    
