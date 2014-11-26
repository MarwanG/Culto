/*
 * CULTO    : Concurrent UPMC LapTop Orchestra
 *            Exercice PPC - 2014
 *
 * Exercice 2. Création d'un objet PureData de compteur borné donnant des listes de "pouet" en sortie
 *
 * author   : Philippe Esling
 *            <esling@ircam.fr>
 *
 */

#ifndef _MULTIPOUET_H_
# define _MULTIPOUET_H_

# include "m_pd.h"

static t_class  *multipouet_class;
typedef struct  _multipouet
{
    t_object    x_obj;
    t_int       i_count;
    t_float     step;
    t_int       i_min, i_max;
    t_outlet    *p_out, *b_out;
}
t_multipouet;

/*
 * Q.4 - Comportement de l'objet en cas de messages reset, set ou bound
 */
void            multipouet_reset(t_multipouet *x);
void            multipouet_set(t_multipouet *x, t_floatarg f);
void            multipouet_bound(t_multipouet *x, t_floatarg min, t_floatarg max);

/*
 * Q.3 - Comportement en cas de message bang
 */
void            multipouet_bang(t_multipouet *x);

/*
 * Q.2 - Création d'un nouvel objet multipouet
 */
void            *multipouet_new(t_symbol *s, int argc, t_atom *argv);

/*
 * Q.1 - Chargement en mémoire des objets de type multipouet
 */
void            multipouet_setup(void);

#endif