import numpy as np
import matplotlib.pyplot as plt

zi = complex(0, 1)

def windowFFT(tgrid, \
              field, \
              nf, lf, \
              tt = ' ', \
              show_plots = False, \
              plot_format = 'display'):
    nf = int(nf)
    lf = float(lf)
    tgrid = np.array(tgrid)
    tgrid_n = (tgrid - tgrid[0]) / (tgrid[-1] - tgrid[0])
    window = np.cos(np.pi * tgrid_n - np.pi / 2.)
    if 1 == 0:
        plt.plot(tgrid, abs(field), label = 'abs')
        plt.plot(tgrid, np.real(field), '+-', label = 'real')
        plt.plot(tgrid, np.imag(field), '+-', label = 'imag')
        plt.plot(tgrid, window)
        plt.title(tt)
        plt.legend()
        plt.show()

    fgrid = np.linspace(-lf,lf,nf,endpoint = False)
    field_f = np.empty(0, dtype = 'complex128')
    for f in fgrid:
        numerator = 0.
        denominator = 0.
        for i in range(len(field) - 1):
            numerator += 0.5*(field[i] * window[i] * \
                         np.exp(zi*f*tgrid[i]) + \
                         field[i + 1] * window[i + 1] * \
                         np.exp(zi*f*tgrid[i + 1])) * \
                          (tgrid[i + 1] - tgrid[i])
            denominator += 0.5 * (window[i] + window[i + 1]) * \
                          (tgrid[i + 1] - tgrid[i])
        field_f = np.append(field_f, numerator / denominator)
    if show_plots:
        plt.figure()
        plt.plot(fgrid, abs(field_f), label = 'abs')
        plt.plot(fgrid, np.real(field_f), '+-', label = 'real')
        plt.plot(fgrid, np.imag(field_f), '+-', label = 'imag')
        plt.xlabel('f (cref/Lref)')
        plt.legend()
        plt.title(tt)
        if plot_format == 'display':
            plt.show()
        elif plot_format == 'ps':
            filename = 'f_' + tt + '.ps'
            fig=plt.gcf()
            fig.savefig(filename, format = 'ps', bbox_inches = 'tight')
    return fgrid, field_f
