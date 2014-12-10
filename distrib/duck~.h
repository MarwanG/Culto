/*
 * CULTO    : Concurrent UPMC LapTop Orchestra
 *            Exercice PPC - 2014
 *
 * Exercice 3. Création d'un objet PureData de ducking (modification du volume d'un signal par le volume d'un autre signal).
 *
 * author   : Philippe Esling
 *            <esling@ircam.fr>
 *
 */

#ifndef _DUCK_TILDE_H_
# define _DUCK_TILDE_H_

# include "m_pd.h"

static t_class      *duck_tilde_class;
typedef struct      _duck_tilde
{
    t_object        x_obj;
    t_sample		f_pan;
    t_sample		f;
    t_inlet 		*x_in2;
    t_outlet		*x_out;
}
t_duck_tilde;

/*
 * Q.5 - Fonction centrale effectuant le calcul
 */
t_int           *duck_tilde_perform(t_int *w);

/*
 * Q.4 - Ajout de l'objet duck~ à l'arbre de traitement DSP
 */
void            duck_tilde_dsp(t_duck_tilde *x, t_signal **sp);

/*
 * Q.3 - Libération de la mémoire de l'objet duck~
 */
void            duck_tilde_free(t_duck_tilde *x);

/*
 * Q.2 - Création d'un nouvel objet duck~
 */
void            *duck_tilde_new(void);

/*
 * Q.1 - Chargement en mémoire des objets de type duck~
 */
void            duck_tilde_setup(void);

#endif