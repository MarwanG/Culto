#include "m_pd.h"
#include "horloge.h"
#include <time.h>


void * horloge_new(void){
    t_horloge *h;
    h =(t_horloge*) pd_new(horloge_class);
    h->h_out = outlet_new(&h->x_obj, &s_symbol);
    return (void*)h;


}


void horloge_bang(t_horloge *x){
    time_t rawtime;
    struct tm * timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    outlet_symbol(x->h_out, gensym(asctime(timeinfo)));

}

void horloge_setup(void){
    horloge_class = class_new(gensym("horloge"), 
                    (t_newmethod) horloge_new,
                    0, sizeof(t_horloge), CLASS_DEFAULT, 0);

    class_addbang(horloge_class, horloge_bang);

}


