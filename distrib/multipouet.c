#include "multipouet.h"
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

void * multipouet_new(t_symbol *s, int argc, t_atom *argv){
    t_multipouet *m;
    m =(t_multipouet*) pd_new(multipouet_class);
    m->p_out = outlet_new(&m->x_obj, &s_symbol);
    m->step = 1;
    t_float f1=0, f2=0;
    
    switch(argc){
        default:
        case 3:
            m->step=atom_getfloat(argv+2);
        case 2:
            f2=atom_getfloat(argv+1);
        case 1:
            f1=atom_getfloat(argv);
            break;
        case 0:
            break;
    }

    if(argc<2)
        f2=f1;
    if(f2 < f1){
        m->i_min = f2;
        m->i_max = f1;
    }else{
        m->i_min = f1;
        m->i_max = f2;
    }
    m->i_count = m->i_min; 

    inlet_new(&m->x_obj, &m->x_obj.ob_pd,
            gensym("list"), gensym("bound"));
    floatinlet_new(&m->x_obj, &m->step);

    m->p_out = outlet_new(&m->x_obj, &s_float);
    m->b_out = outlet_new(&m->x_obj, &s_bang);


    return (void*)m;


}


void multipouet_bang(t_multipouet *x){
    while(x->i_count < x->i_min)
        x->i_count = x->i_count+1;

    if(x->i_count > x->i_max)
        x->i_count = x->i_min;
    char *print_pouet = NULL;
    print_pouet = malloc(1+strlen("pouet ")*(x->i_count)*sizeof(char));
    if(print_pouet == NULL){
        perror("malloc");
        exit(EXIT_FAILURE);
    }
        
    int i;
    *print_pouet = '\0';
    for(i=0; i<x->i_count; i++){
        strcat(print_pouet,"pouet ");
    }
    strcat(print_pouet,"\0");
    
    outlet_float(x->p_out,x->i_count);
    outlet_symbol(x->p_out, gensym((const char*)print_pouet));
    x->i_count = x->i_count+x->step;

    printf("%dok?%s\n",x->i_count,print_pouet);
    free(print_pouet);

}

void multipouet_set(t_multipouet *x, t_floatarg f){
    x->i_count = f;
    printf("%f\n",f);
}

void multipouet_bound(t_multipouet *x, t_floatarg min, t_floatarg max){
    
    x->i_min=min;
    x->i_max=max;
    multipouet_reset(x);
}

void multipouet_reset(t_multipouet *x){
    x->i_count = x->i_min;
}

void multipouet_setup(void){
    multipouet_class = class_new(gensym("multipouet"), 
                    (t_newmethod) multipouet_new,
                    0, sizeof(t_multipouet),
                    CLASS_DEFAULT, A_GIMME,0);

    class_addbang(multipouet_class, multipouet_bang);
    class_addmethod(multipouet_class, 
        (t_method)multipouet_set, 
        gensym("set"),A_DEFFLOAT,0);
    class_addmethod(multipouet_class, 
        (t_method)multipouet_reset, 
        gensym("reset"),A_DEFFLOAT,0);
    class_addmethod(multipouet_class, 
        (t_method)multipouet_bound, 
        gensym("bound"),A_DEFFLOAT, A_DEFFLOAT,0);

}


