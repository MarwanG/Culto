/*
 * CULTO    : Concurrent UPMC LapTop Orchestra
 *            Exercice PPC - 2014
 *
 * Exercice 1. Création d'un objet PureData donnant l'heure
 *
 * author   : Philippe Esling
 *            <esling@ircam.fr>
 *
 */

#ifndef _HORLOGE_H_
# define _HORLOGE_H_

# include "m_pd.h"

static t_class 	*horloge_class;
typedef struct 	_horloge
{
    t_object	x_obj;
    t_outlet	*h_out;
}
t_horloge;

/*
 * Q.4 - Comportement de l'objet en cas de message bang
 */
void            horloge_bang(t_horloge *x);

/*
 * Q.3 - Création d'un nouvel objet horloge
 */
void            *horloge_new(void);

/*
 * Q.2 - Chargement en mémoire des objets de type horloge
 */
void            horloge_setup(void);

#endif