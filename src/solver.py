import numpy as np

def build_up_b(dx, dy,
               u, v,
               rho,
               dt
               ):
    """
    Calcola il termine noto :math:`b` dell'equazione di Poisson per la pressione.

    Questo termine garantisce l'incomprimibilità del campo di velocità, forzando
    la divergenza a zero al passo temporale successivo. Per ulteriori dettagli sulla
    derivazione fisica e matematica del termine :math:`b`, consultare la sezione
    :ref:`incompressible_fluid`.

    Parameters
    ----------
    dx : float
        Spaziatura della griglia nella direzione x.
    dy : float
        Spaziatura della griglia nella direzione y.
    rho : float
        Densità del fluido.
    dt : float

    Returns
    -------

    """
    b = np.zeros_like(u)

    b[1:-1, 1:-1] = (rho * (1 / dt * 
                    ((u[1:-1, 2:] - u[1:-1, 0:-2]) / 
                     (2 * dx) + (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -
                    ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 -
                      2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) *
                           (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx))-
                          ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2))

    return b