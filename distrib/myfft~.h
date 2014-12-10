/*
 * CULTO    : Concurrent UPMC LapTop Orchestra
 *            Exercice PPC - 2014
 *
 * Exercice 4. Calculer la FFT par soi-meme et comparer ses résultats.
 *
 * author   : Philippe Esling
 *            <esling@ircam.fr>
 *
 */

#ifndef _myfft_TILDE_H_
# define _myfft_TILDE_H_

# include "m_pd.h"

#ifndef PI
#define PI 3.141592653589793115997963468544185161590576171875
#endif

#ifndef TWOPI
#define TWOPI 6.28318530717958623199592693708837032318115234375
#endif

static t_class      *myfft_tilde_class;
typedef struct      _myfft_tilde
{
    t_object        x_obj;
    t_sample		f_pan;
    t_sample		f;
    t_outlet		*x_out;
    int             *bitshuffle;
    float           *weighting;
    float           *window;
}
t_myfft_tilde;

/*
 * Q.5 - Fonction centrale effectuant le calcul
 */
t_int           *myfft_tilde_perform(t_int *w);

/*
 * Q.4 - Ajout de l'objet myfft~ à l'arbre de traitement DSP
 */
void            myfft_tilde_dsp(t_myfft_tilde *x, t_signal **sp);

/*
 * Q.3 - Libération de la mémoire de l'objet myfft~
 */
void            myfft_tilde_free(t_myfft_tilde *x);

/*
 * Q.2 - Création d'un nouvel objet myfft~
 */
void            *myfft_tilde_new(void);

/*
 * Q.1 - Chargement en mémoire des objets de type myfft~
 */
void            myfft_tilde_setup(void);

#endif