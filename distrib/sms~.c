#include "sms~.h"
#include "myfft~.h"
#include "myfft_fft.h"
#include <stdlib.h>
#include <math.h>
#define BUFFER_LEN 4096

void fen_hamming(t_sample* buf, int size){
    int i;
    for(i=0;i<size;i++){
      if(buf[i] < size || buf[i] < 0){
	    buf[i] = 0.54 - (0.46 * cos(2 * PI *(buf[i]/size)));
      }else{
	    buf[i] = 0;
      }
    }
}


t_int *sms_tilde_perform (t_int *w){
    printf("Starting Perform...\n");
    t_sms_tilde* x = (t_sms_tilde*) w[1];
    t_sample* in_srcHarm = (t_sample*) w[2];
    t_sample* in_srcMod = (t_sample*) w[3];
    t_sample* out = (t_sample*) w[4];
    t_sample* bufferHarm = x->bufferHarm;
    t_sample* bufferMod = x->bufferMod;
    int vecSize = (int) w[5];
    int i=0, inGain=0, curGain=0, outGain=0, finalGain=0;
    printf("bp : %d\n",x->bypass);
    printf("ano : %d\n",x->autonorm);
    printf("seuil : %f\n",x->seuil);
    
    /* Buffer gestion */
    for (i = 0; i < vecSize; i++) { // Filling buffers with incoming data
        if (x->cptBufHarm < x->buffer_size && x->cptBufMod < x->buffer_size) {
            bufferHarm[x->cptBufHarm] = in_srcHarm[i];
            bufferMod[x->cptBufMod] = in_srcMod[i];
            x->cptBufHarm++;
            x->cptBufMod++;
        }
    }
    if (x->bypass == 1){ // if bypass
        for (i = 0; i < vecSize; i++) {
            out[i] = in_srcHarm[i];
        }
        int j=0;
        for (i = vecSize; i < x->cptBufHarm; i++) { 
        // removing the $vecSize first data of the buffer
            bufferHarm[j] = bufferHarm[i];
            j++;
        }
        x->cptBufHarm -= vecSize;
        return w+6;
    }
    printf("cptO : %d, vecS : %d\n",x->cptOut,vecSize);
    if (x->cptOut < vecSize){
    //If our output buffer is already ready to out, skip the process
    //The data is stored in other buffers anyway
        printf("Step 1\n");
        /*
        t_sample* dupHarm = (t_sample*) malloc((2+vecSize)*sizeof(t_sample));
        t_sample* dupMod = (t_sample*) malloc((2+vecSize)*sizeof(t_sample));
        for (i = 0; i < vecSize; i++) {
            dupHarm[i] = in_srcHarm[i];
            dupMod[i] = in_srcMod[i];
        }*/
        printf("Step 2\n");
        printf("Step 3\n");
        int minCpt=0, maxCpt=0; // getting the size of the smaller buffer
        if (x->cptBufHarm > x->cptBufMod){
            minCpt = x->cptBufMod;
            maxCpt = x->cptBufHarm;
        }else {
            minCpt = x->cptBufHarm;
            maxCpt = x->cptBufMod;
        }
            
        fen_hamming(bufferHarm, minCpt);
        fen_hamming(bufferMod, minCpt);

        int* bitshuffle=NULL;
        bitshuffle = (int*)malloc((2*minCpt)*sizeof(int));
        float* weightning=NULL;
        weightning = (float*)malloc((2*minCpt)*sizeof(float));
        init_rdft(minCpt, bitshuffle, weightning);
        printf("Step 4\n");
        rdft(minCpt, 1, bufferHarm, bitshuffle, weightning);
        rdft(minCpt, 1, bufferMod, bitshuffle, weightning);
        printf("Step 5\n");


        if (x->autonorm == 1){ // if autonorm
            inGain=0;
            for (i = 0; i < minCpt-1; i+=2) {
                inGain += fabsf(bufferHarm[i]-bufferHarm[i+1]);
            }
        } 
        float* tmp1=NULL; 
        float* tmp2=NULL;
        tmp1 = (float*)malloc(x->buffer_size*sizeof(float));
        if(tmp1 == NULL){
            perror("malloc");
            exit(EXIT_FAILURE);
        }
        tmp2 = (float*)malloc(x->buffer_size*sizeof(float));
        if(tmp2 == NULL){
            perror("malloc");
            exit(EXIT_FAILURE);
        }
        
        for (i = 0; i < minCpt-1; i+=2) {
            curGain = fabsf(bufferMod[i]-bufferMod[i+1]);
            if (curGain > x->seuil) {
                tmp1[i] = fabsf(bufferHarm[i]-bufferHarm[i+1]);
            }
            tmp1[i+1] = - atan2(bufferHarm[i+1], bufferHarm[i]);
            tmp2[i] = tmp1[i] * cos(tmp1[i+1]);
            tmp2[i+1] = -tmp1[i] * sin(tmp1[i+1]);
        }
        
        if (x->autonorm == 1) {
            outGain = 0;
            for (i = 0; i < minCpt-1; i+=2) {
                outGain += fabsf(tmp2[i]-tmp2[i+1]);
            }
            finalGain = inGain/outGain;
        }else {
            finalGain = 1;
        }
        rdft(minCpt, -1, tmp2, bitshuffle, weightning);

        for (i = 0; i < minCpt; i++) { // filling the output buffer
            if(x->cptOut < x->buffer_size){
                tmp2[i] = tmp2[i]/finalGain;
                x->bufferOut[x->cptOut] = tmp2[i];
                x->cptOut++;
                printf("test\n");
            }
        }
        // We empty the buffers, data is already processed and stored in bufferOut
        int j=0;
        for (i = minCpt; i < x->cptBufHarm; i++) {
            bufferHarm[j] = bufferHarm[i];
            j++;
        }
        x->cptBufHarm -= minCpt;
        j=0;
        for (i = minCpt; i< x->cptBufMod; i++){
            bufferMod[j] = bufferMod[i];
            j++;
        }
        x->cptBufMod -= minCpt;

        free(tmp1);
        free(tmp2);
        free(bitshuffle);
        free(weightning);
    }
    // Filling the real output vector
    for (i = 0; i < vecSize; i++) {
        out[i] = x->bufferOut[i];
    }
    int j=0;
    for (i = vecSize; i < x->cptOut; i++) {
        x->bufferOut[j] = x->bufferOut[i];
        j++;
    }
    printf("coucou %d\n",x->cptOut);
    x->cptOut -= vecSize;

    printf("coucou %d\n",x->cptOut);
    return w+6;
}

void sms_tilde_dsp(t_sms_tilde *x,t_signal **sp){
    dsp_add(sms_tilde_perform, 5, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[0]->s_n); 

}

void sms_tilde_free(t_sms_tilde *x){
	inlet_free(x->in1);
	inlet_free(x->in2);
	inlet_free(x->in3);
	inlet_free(x->in4);
	outlet_free(x->out);
	free(x->bufferHarm);
	free(x->bufferMod);
}

void sms_tilde_msg(t_sms_tilde *x, t_floatarg auton, t_floatarg bypa){
    x->autonorm = auton;
    x->bypass = bypa;
}

void *sms_tilde_new(t_symbol *s, int argc, t_atom *argv){
	t_sms_tilde *m;
	m = (t_sms_tilde*)pd_new(sms_tilde_class);
	m->in1 = inlet_new(&m->x_obj,&m->x_obj.ob_pd,&s_signal,&s_signal);
	m->in2 = inlet_new(&m->x_obj,&m->x_obj.ob_pd,&s_signal,&s_signal);
 	m->in3 = floatinlet_new(&m->x_obj,&m->seuil);
    inlet_new(&m->x_obj, &m->x_obj.ob_pd,
            gensym("list"), gensym("msg"));
 	m->out = outlet_new(&m->x_obj,&s_signal);
    m->bufferOut = (t_sample*)malloc(BUFFER_LEN*sizeof(t_sample));

 	//if(argc < 2){
 		m->bufferHarm = (t_sample*)malloc(BUFFER_LEN * sizeof (t_sample));
 		m->bufferMod = (t_sample*)malloc(BUFFER_LEN * sizeof (t_sample));
        m->buffer_size = BUFFER_LEN;
 	/*}else{
 		m->buffer = (t_sample*)malloc(atom_getint(argv) * sizeof (t_sample));
        m->buffer_size = atom_getint(argv);
 	}*/
    m->cptBufMod = 0;
    m->cptBufHarm = 0;
    m->cptOut = 0;
        printf("%d\n",m->buffer_size);

 	return (void*) m;
}

void sms_tilde_setup(void){
	sms_tilde_class = class_new(gensym("sms~"),(t_newmethod)sms_tilde_new
		,(t_method)sms_tilde_free,sizeof(t_sms_tilde),
		CLASS_DEFAULT,0);
	class_addmethod(sms_tilde_class, (t_method)sms_tilde_dsp,gensym("dsp"),0);
    class_addmethod(sms_tilde_class, 
        (t_method)sms_tilde_msg, 
        gensym("msg"),A_DEFFLOAT, A_DEFFLOAT,0);
   	CLASS_MAINSIGNALIN(sms_tilde_class, t_sms_tilde, f);  

}
void            init_rdft(int n, int *ip, float *w)
{
  int           nw,nc;
  printf("init rdft\n");
  nw = n >> 2;
  printf("1\n");
  makewt(nw, ip, w);
  printf("2\n");
  nc = n >> 2;
  printf("3\n");
  makect(nc, ip, w + nw);
  printf("4\n");
  return;
}

/*
 *  @brief      Computes the FFT based on discrete (Dft) and real (Rdft) numbers
 *
 *  @param      n       : Size of the vector to transform
 *  @param      isgn    : Flag for signed (1) or unsigned (0) data
 *  @param      a       : Buffer of data to transform
 *  @param      ip      : Bitshuffling vector of size n*2 (to be inited with init_rdft)
 *  @param      ip      : Weighting vector of size n*2 (to be inited with init_rdft)
 *
 */
void            rdft(int n, int isgn, float *a, int *ip, float *w)
{

  int           j, nw, nc;
  float         xi;
  void          bitrv2(int n, int *ip, float *a),
                cftsub(int n, float *a, float *w),
                rftsub(int n, float *a, int nc, float *c);
  nw = ip[0];
  nc = ip[1];
  if (isgn < 0)
  {
    a[1] = 0.5 * (a[1] - a[0]);
    a[0] += a[1];
    for (j = 3; j <= n - 1; j += 2)
      a[j] = -a[j];

    if (n > 4) {
      rftsub(n, a, nc, w + nw);
      bitrv2(n, ip + 2, a);
    }

    cftsub(n, a, w);

    for (j = 1; j <= n - 1; j += 2) {
      a[j] = -a[j];
    }
  }

  else {

    if (n > 4) {
      bitrv2(n, ip + 2, a);
    }

    cftsub(n, a, w);

    if (n > 4) {
      rftsub(n, a, nc, w + nw);
    }

    xi = a[0] - a[1];
    a[0] += a[1];
    a[1] = xi;
  }
}


void            bitrv2(int n, int *ip, float *a)
{
  int           j, j1, k, k1, l, m, m2;
  float         xr, xi;
    
  ip[0] = 0;
  l = n;
  m = 1;
  while ((m << 2) < l)
  {
    l >>= 1;
    for (j = 0; j <= m - 1; j++) {
      ip[m + j] = ip[j] + l;
    }
    m <<= 1;
  }
  if ((m << 2) > l) {

    for (k = 1; k <= m - 1; k++) {

      for (j = 0; j <= k - 1; j++) {
	j1 = (j << 1) + ip[k];
	k1 = (k << 1) + ip[j];
	xr = a[j1];
	xi = a[j1 + 1];
	a[j1] = a[k1];
	a[j1 + 1] = a[k1 + 1];
	a[k1] = xr;
	a[k1 + 1] = xi;
      }
    }
  }

  else {
    m2 = m << 1;

    for (k = 1; k <= m - 1; k++) {

      for (j = 0; j <= k - 1; j++) {
	j1 = (j << 1) + ip[k];
	k1 = (k << 1) + ip[j];
	xr = a[j1];
	xi = a[j1 + 1];
	a[j1] = a[k1];
	a[j1 + 1] = a[k1 + 1];
	a[k1] = xr;
	a[k1 + 1] = xi;
	j1 += m2;
	k1 += m2;
	xr = a[j1];
	xi = a[j1 + 1];
	a[j1] = a[k1];
	a[j1 + 1] = a[k1 + 1];
	a[k1] = xr;
	a[k1 + 1] = xi;
      }
    }
  }
}


void        cftsub(int n, float *a, float *w)
{
  int       j, j1, j2, j3, k, k1, ks, l, m;
  float     wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
  float     x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    
  l = 2;

  while ((l << 1) < n) {
    m = l << 2;

    for (j = 0; j <= l - 2; j += 2) {
      j1 = j + l;
      j2 = j1 + l;
      j3 = j2 + l;
      x0r = a[j] + a[j1];
      x0i = a[j + 1] + a[j1 + 1];
      x1r = a[j] - a[j1];
      x1i = a[j + 1] - a[j1 + 1];
      x2r = a[j2] + a[j3];
      x2i = a[j2 + 1] + a[j3 + 1];
      x3r = a[j2] - a[j3];
      x3i = a[j2 + 1] - a[j3 + 1];
      a[j] = x0r + x2r;
      a[j + 1] = x0i + x2i;
      a[j2] = x0r - x2r;
      a[j2 + 1] = x0i - x2i;
      a[j1] = x1r - x3i;
      a[j1 + 1] = x1i + x3r;
      a[j3] = x1r + x3i;
      a[j3 + 1] = x1i - x3r;
    }

    if (m < n) {
      wk1r = w[2];

      for (j = m; j <= l + m - 2; j += 2) {
	j1 = j + l;
	j2 = j1 + l;
	j3 = j2 + l;
	x0r = a[j] + a[j1];
	x0i = a[j + 1] + a[j1 + 1];
	x1r = a[j] - a[j1];
	x1i = a[j + 1] - a[j1 + 1];
	x2r = a[j2] + a[j3];
	x2i = a[j2 + 1] + a[j3 + 1];
	x3r = a[j2] - a[j3];
	x3i = a[j2 + 1] - a[j3 + 1];
	a[j] = x0r + x2r;
	a[j + 1] = x0i + x2i;
	a[j2] = x2i - x0i;
	a[j2 + 1] = x0r - x2r;
	x0r = x1r - x3i;
	x0i = x1i + x3r;
	a[j1] = wk1r * (x0r - x0i);
	a[j1 + 1] = wk1r * (x0r + x0i);
	x0r = x3i + x1r;
	x0i = x3r - x1i;
	a[j3] = wk1r * (x0i - x0r);
	a[j3 + 1] = wk1r * (x0i + x0r);
      }

      k1 = 1;
      ks = -1;

      for (k = (m << 1); k <= n - m; k += m) {
	k1++;
	ks = -ks;
	wk1r = w[k1 << 1];
	wk1i = w[(k1 << 1) + 1];
	wk2r = ks * w[k1];
	wk2i = w[k1 + ks];
	wk3r = wk1r - 2 * wk2i * wk1i;
	wk3i = 2 * wk2i * wk1r - wk1i;

	for (j = k; j <= l + k - 2; j += 2) {
	  j1 = j + l;
	  j2 = j1 + l;
	  j3 = j2 + l;
	  x0r = a[j] + a[j1];
	  x0i = a[j + 1] + a[j1 + 1];
	  x1r = a[j] - a[j1];
	  x1i = a[j + 1] - a[j1 + 1];
	  x2r = a[j2] + a[j3];
	  x2i = a[j2 + 1] + a[j3 + 1];
	  x3r = a[j2] - a[j3];
	  x3i = a[j2 + 1] - a[j3 + 1];
	  a[j] = x0r + x2r;
	  a[j + 1] = x0i + x2i;
	  x0r -= x2r;
	  x0i -= x2i;
	  a[j2] = wk2r * x0r - wk2i * x0i;
	  a[j2 + 1] = wk2r * x0i + wk2i * x0r;
	  x0r = x1r - x3i;
	  x0i = x1i + x3r;
	  a[j1] = wk1r * x0r - wk1i * x0i;
	  a[j1 + 1] = wk1r * x0i + wk1i * x0r;
	  x0r = x1r + x3i;
	  x0i = x1i - x3r;
	  a[j3] = wk3r * x0r - wk3i * x0i;
	  a[j3 + 1] = wk3r * x0i + wk3i * x0r;
	}
      }
    }

    l = m;
  }

  if (l < n) {

    for (j = 0; j <= l - 2; j += 2) {
      j1 = j + l;
      x0r = a[j] - a[j1];
      x0i = a[j + 1] - a[j1 + 1];
      a[j] += a[j1];
      a[j + 1] += a[j1 + 1];
      a[j1] = x0r;
      a[j1 + 1] = x0i;
    }
  }
}


void        rftsub(int n, float *a, int nc, float *c)
{
  int       j, k, kk, ks;
  float     wkr, wki, xr, xi, yr, yi;
    
  ks = (nc << 2) / n;
  kk = 0;

  for (k = (n >> 1) - 2; k >= 2; k -= 2) {
    j = n - k;
    kk += ks;
    wkr = 0.5 - c[kk];
    wki = c[nc - kk];
    xr = a[k] - a[j];
    xi = a[k + 1] + a[j + 1];
    yr = wkr * xr - wki * xi;
    yi = wkr * xi + wki * xr;
    a[k] -= yr;
    a[k + 1] -= yi;
    a[j] += yr;
    a[j + 1] -= yi;
  }
}


void        makewt(int nw, int *ip, float *w)
{
    int     nwh, j;
    float   delta, x, y;
    
    ip[0] = nw;
    ip[1] = 1;
    if (nw > 2) {
        nwh = nw >> 1;
        delta = atan(1.0) / nwh;
        w[0] = 1;
        w[1] = 0;
        w[nwh] = cos(delta * nwh);
        w[nwh + 1] = w[nwh];
        for (j = 2; j <= nwh - 2; j += 2) {
            x = cos(delta * j);
            y = sin(delta * j);
            w[j] = x;
            w[j + 1] = y;
            w[nw - j] = y;
            w[nw - j + 1] = x;
        }
        bitrv2(nw, ip + 2, w);
    }
}


void        makect(int nc, int *ip, float *c)
{
    int     nch, j;
    float   delta;
    
    ip[1] = nc;
    if (nc > 1) {
        nch = nc >> 1;
        delta = atan(1.0) / nch;
        c[0] = 0.5;
        c[nch] = 0.5 * cos(delta * nch);
        for (j = 1; j <= nch - 1; j++) {
            c[j] = 0.5 * cos(delta * j);
            c[nc - j] = 0.5 * sin(delta * j);
        }
    }
}

