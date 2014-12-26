#ifndef SMS_DEF
#define SMS_DEF

#include "m_pd.h"

static t_class *sms_tilde_class;

typedef struct _sms_tilde
{
	 t_object   x_obj;
	 t_inlet    *in1;
	 t_inlet    *in2;
	 t_inlet    *in3;
	 t_inlet    *in4;
	 int        buffer_size;
	 t_float    seuil;
	 t_float    info;
	 t_outlet   *out;
	 t_sample   *buffer;
	 t_sample   f;
}t_sms_tilde;



t_int *sms_tilde_perform (t_int *w);

void sms_tilde_dsp(t_sms_tilde *x,t_signal **sp);

void sms_tilde_free(t_sms_tilde *x);

void *sms_tilde_new(t_symbol *s, int argc, t_atom *argv);

void sms_tilde_setup(void);

#endif
