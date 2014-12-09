#include "m_pd.h"

static t_class *duck_tilde_class;
typedef struct _duck_tilde
{
  t_object x_obj;
  t_sample f_pan;
  t_sample f;
  t_inlet *x_in2;
  t_outlet *x_out;
}
t_duck_tilde;

t_int *duck_tilde_perform (t_int *w);

void duck_tilde_dsp(t_duck_tilde *x,t_signal **sp);

void duck_tilde_free(t_duck_tilde *x);

void *duck_tilde_new(void);

void duck_tilde_setup(void);
