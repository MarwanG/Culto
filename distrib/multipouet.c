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
  if(x->i_count > x->i_max){
    multipouet_reset(x);
    outlet_bang(x->b_out);
  }
   t_atom list[x->i_count];
   for(i=0;i<x->i_count;i++){
     SETSYMBOL(list,gensym("pouet"));
   }
   outlet_list(x->p_out,gensym("list"),x->i_count,list);
}

/*
 * Q.2 - Création d'un nouvel objet multipouet
 */
void            *multipouet_new(t_symbol *s, int argc, t_atom *argv) /* ordre arguments : min, max, step */
{
    t_multipouet   *m;
    m = (t_multipouet *)pd_new(multipouet_class);
    if(argc ==3){
      m->i_min = atom_getintarg(0,argc,argv);
      m->i_max =  atom_getintarg(1,argc,argv);
      m->step =  atom_getintarg(2,argc,argv);
    }
    else if (argc == 2){
      m->i_min =  atom_getintarg(0,argc,argv);
      m->i_max =  atom_getintarg(1,argc,argv);
      m->step = 1;
    }
    else if (argc == 1){
      m->i_min =  atom_getintarg(0,argc,argv);
      m->i_max = 10;
      m->step = 1;
    }
    else{
      m->step = 1; 
      m->i_min = 0;
      m->i_max = 10;
    }
    m->i_count = 0;
    m->p_out = outlet_new(&m->x_obj,&s_symbol);
    m->b_out = outlet_new(&m->x_obj,&s_bang);
    floatinlet_new(m,&m->step); 
    inlet_new(m,&m->x_obj,gensym("list"),gensym("bound")); 
    return (void *)m;
}

/*
 * Q.1 - Chargement en mémoire des objets de type multipouet
 */
void multipouet_setup(void){
  multipouet_class = class_new(gensym("multipouet"),
                              (t_newmethod)multipouet_new,
                              0, sizeof(t_multipouet),
         CLASS_PATCHABLE, A_GIMME, 0 );
  class_addbang(multipouet_class,multipouet_bang);
	class_addmethod(multipouet_class, multipouet_reset,gensym("reset"),0);
  class_addmethod(multipouet_class, multipouet_set,gensym("set"),A_DEFFLOAT,0);
  class_addmethod(multipouet_class, multipouet_bound,gensym("bound"),A_DEFFLOAT,A_DEFFLOAT,0);		 
}
