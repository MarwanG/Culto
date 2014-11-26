#include "m_pd.h"
#include "multipouet.h"





void multipouet_reset(t_multipouet *x){
  x->i_count = x->i_min;
}
void multipouet_set(t_multipouet *x, t_floatarg f){
  x->i_count = f;
}
void multipouet_bound(t_multipouet *x, t_floatarg min, t_floatarg max){
  x->i_min = min;
  x->i_max = max;
}

/*
 * Q.3 - Comportement en cas de message bang
 */
void multipouet_bang(t_multipouet *x){
  int i;
  x->i_count = x->i_count + 1 ;
  if(x->i_count > i_max){
    multipouet_reset(x);
    outlet_bang(x->b_out);
  }
   const char *list[i_count];
   for(i=0;i<x->i_count;i++){
     list[i]="pouet"; 
   }
   outlet_list(x->p_out,gensysm("list"),x->i_count,list);
}

/*
 * Q.2 - Création d'un nouvel objet multipouet
 */
void *multipouet_new(t_symbol *s, int argc, t_atom *argv){
  t_multipouet * mp;
  mp = (t_multipouet*) pd_new (multipouet_class);
  mp->i_min =atom_getint(argv[0]);
  mp->i_max = atom_getint(argv[1]);
  mp->step = atom_getfloat(argv[2]);
  mp->i_count = mp->i_min;
  mp->p_out=outlet_new(&mp->x_obj,&s_list);
  mp->b_out=outlet_new(&mp->x_obj,&s_bang);

  /* pas sur pour l'inlet je le comprend pas pour l'instant */
  
  return (void*)mp;
}

/*
 * Q.1 - Chargement en mémoire des objets de type multipouet
 */
void multipouet_setup(void){
  multipouet_class = class_new(gensym("multipouet"),
			       (t_newmethod) multipouet_new,
			       0,sizeof(t_multipouet),CLASS_DEFAULT,A_GIMME);
  class_add(multipouet_class,multipouet_bang);
			   
}
