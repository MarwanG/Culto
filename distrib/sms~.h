#include "m_pd.h"

static t_class *sms_tilde_class;

typedef struct _sms_tilde
{
	 t_object x_obj;
	 t_inlet *x_in1;
	 t_inlet *x_in2;
	 t_inlet *x_in3;
	 t_inlet *x_in4;
	 int buffer_size;
	 t_float seuil;
	 t_gpointer info;
	 t_outlet *x_out;
	 t_sample *buffer;
}t_sms_tilde;


t_int *sms_tilde_perform (t_int *w);

void sms_tilde_dsp(t_duck_tilde *x,t_signal **sp);

void sms_tilde_free(t_duck_tilde *x);

void *sms_tilde_new(void);

void sms_tilde_setup(void);