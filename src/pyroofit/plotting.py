# -*- coding: utf-8 -*-
""" Plot function for the PDF class

Tools to make nice plots from RooFit pdfs.
The function fast_plot is used by the PDF class
to make default plots.

Todo:
    * Plotter class containing the RooFit frame
    * Equal length ticks
    * Provide matplolib functionality


"""
from __future__ import print_function
from .utilities import  ClassLoggingMixin
from .data import th12uarray, roo2hist
import ROOT

import numpy as np
import root_numpy
from uncertainties import ufloat
from uncertainties import unumpy as unp

import matplotlib.pyplot as plt
#from matplotlib.ticker import NullFormatter
import matplotlib


class Plotter(ClassLoggingMixin):
    """ Experimental Plotter class

    This function serves the purpose to create a RooFit frame without the need to interface
    RooFit.

    Todo:
        * Adding plot pdf functionality

    """
    def __init__(self, pdf, observable=None, nbins=20):
        super(Plotter, self).__init__()
        self.pdf = pdf
        self.observable = observable if observable is not None else pdf.get_observable()
        self.frame = None
        self.create_frame()
        self.nbins = nbins

    def create_frame(self, title="Fit"):
        self.frame = self.pdf.get_observable().frame(ROOT.RooFit.Title(title), ROOT.RooFit.Bins(self.nbins))


def get_optimal_bin_size(n, round=True):
    """Helper function to calculate optimal binning

    This function calculates the optimal amount of bins for the number of events n.

    Args:
        n (int): number of events to be binned
        round (bool or int): Round to
    Returns:
        (int): Optimal number of bins

    """

    def roundtobase(n, base=5):
        diff = n % base
        if diff <= base / 2.:
            return n - diff
        else:
            return n - diff + base

    n_opt = int(2 * n**(1/3.0))

    if round:
        base = 5
        if isinstance(round, int):
            base = round
        n_opt = roundtobase(n_opt, base)

    return n_opt


def round_to_1(x):
    from math import log10, floor
    return round(x, -int(floor(log10(abs(x)))))




def py_plot(model, data, observable, filename=None, components=None, title = None,
            nbins = None, round_bins = 5, legend = True,
            legend_kwargs = {}, figure_kwargs = {}, title_kwargs = {},
            xlabel_kwargs = {}, pull_ylabel_kwargs = {}, ylabel_kwargs = {},
            legend_data_name="Data", legend_fit_name="Fit",
            fit_color = 'black', data_color = 'black', show=False):
    """ Generic plot function.

    Args:
        model (RooAbsPDF):
            Fit model to be drawn
        data (RooDataSet):
            Dataset to be plotted
        observable (RooAbsVar):
            Observable to be drawn
        filename (str):
            Name of the output file. Suffix determines file type
        components (list(tuple(PDF, RooRealVar))):
            input list of pdf components to be plotted (automatic for composite pdfs)
        title (str):
            title for plot
        nbins (int):
            Number of bins
        round_bins (int):
            magic to for automatically choosing the bin numbers
        legend (bool):
            bool to enable or disable plotting of the legend
        legend_kwargs (dict):
            keyword arguments passed to the matplotlib.axes.Axes.legend instance
        figure_kwargs (dict):
            keyword arguments passed to the matplotlib.pyplot.figure instance
        title_kwargs (dict):
            keyword arguments passed to the title
        xlabel_kwargs (dict):
            keyword arguments passed to the xlabel
        ylabel_kwargs (dict):
            keyword arguments passed to the ylabel
        pull_ylabel_kwargs (dict):
            keyword arguments passed to the ylabel of the pull
        legend_data_name (str):
            Name of the data part in the fit plot
        legend_fit_name (str):
            name of the total fit in the plot
        fit_color (str):
            color used to plot the fit line
        data_color (str):
            color used to plot the data
        show (bool):
            flag for plt.show() call
    """

    nbins = get_optimal_bin_size(data.numEntries(), round_bins) if nbins is None else nbins

    if isinstance(data, ROOT.RooDataHist):
        nbins = observable.getBins()

    #convert data to useable format for pyplot
    if isinstance(data, ROOT.RooDataSet):
        data  = roo2hist(data, nbins, observable, 'data_roo_hist')
    data_root_hist = data.createHistogram('data_root_hist', observable)
    data_uarray, bins = th12uarray(data_root_hist)

    bins = np.array(bins)
    bin_centers = (bins[:-1] + bins[1:])/2
    bin_widths  = (bins[1:] - bins[:-1])

    #### PLOTTING ####
    # definitions for the axes - #TODO - make this an option
    rect_hist = [0, 0.265, 1, 0.7]
    rect_pull = [0, 0, 1, 0.235]

    # start with a rectangular Figure
    if 'figsize' not in figure_kwargs:
        figure_kwargs['figsize'] = (5,5)
    fig = plt.figure(**figure_kwargs)


    ax_hist = plt.axes(rect_hist)
    ax_pull = plt.axes(rect_pull, sharex=ax_hist)

    # no labels for the hist plot
    plt.setp(ax_hist.get_xticklabels(), visible=False)
    #nullfmt = NullFormatter() # no labels
    #ax_hist.xaxis.set_major_formatter(nullfmt)

    legend_handles = []
    legend_labels = []

    data_plot = ax_hist.errorbar(bin_centers, unp.nominal_values(data_uarray),
                                 yerr=unp.std_devs(data_uarray), xerr= bin_widths / 2,
                                 color=data_color, fmt='.', capthick=0.0, lw=2.0)

    legend_handles.append(data_plot)
    legend_labels.append(legend_data_name)


    npoints_curve = 10000

    if components is not None:
        for component, ni in components:
            component.ni   = ufloat(ni.getVal(), ni.getError())

            if isinstance(component.roo_pdf, ROOT.RooHistPdf):
                #RooHistPDFs have intrinsic binning so need to treated differently
                h = component.roo_pdf.createHistogram(observable.GetName(), npoints_curve)
                hy, hx = th12uarray(h)
                hy     = np.array(unp.nominal_values(hy))
                hy    /= hy.sum()
                hx_centers = (hx[:-1] + hx[1:])/2

                component.plot = ax_hist.hist(hx_centers, weights=hy * component.ni.n, bins=hx,
                                              histtype= 'step', color=component.color, **component.plot_kwargs)

                if component.color == None:
                    component.color = component.plot[2][0].get_ec()

                if component.fill or component.hatch != False:
                    component.face_color = component.color
                    component.edge_color = component.color

                    #alternative for hatch is a string (ex. '/') not True
                    if component.hatch != False and not component.fill:
                        component.face_color = "None"
                    if component.hatch == False and component.fill:
                        component.edge_color = "None"

                    #want the bar shape of the histogram so adjust hx and hy
                    print(len(hx), len(hy))

                    hx_new = []
                    for i in range(len(hx)):
                        if ((i == 0) or (i==len(hx)-1)):
                            hx_new.append(hx[i])
                        else:
                            hx_new.append(hx[i])
                            hx_new.append(hx[i])
                    hx = hx_new

                    hy_new = []
                    for i in range(len(hy)):
                        hy_new.append(hy[i])
                        hy_new.append(hy[i])
                    hy = hy_new

                    hy = np.array(hy)
                    hx = np.array(hx)

                    print(len(hx), len(hy))

                    component.fill_plot = ax_hist.fill_between(hx, np.zeros(len(hy)), hy * component.ni.n,
                    facecolor=component.face_color, edgecolor = component.edge_color,
                    alpha=component.fill_alpha, hatch=component.hatch, **component.fill_kwargs)
                    legend_handles.append( (component.plot[2][0], component.fill_plot) )
                else:
                    legend_handles.append(component.plot[2][0])

                legend_labels.append(component.title)

                continue


            #### all plots but hist
            hx, hy = component.get_curve(observable.GetName(), npoints_curve)
            hy    *= component.ni.n / nbins

            component.plot = ax_hist.plot(hx, hy, color=component.color, **component.plot_kwargs)
            component.fill_plot = None

            if component.color == None:
                component.color = component.plot[0].get_color()

            if component.fill or component.hatch != False:
                component.face_color = component.color
                component.edge_color = component.color

                #alternative for hatch is a string (ex. '/') not True
                if component.hatch != False and not component.fill:
                    component.face_color = "None"
                if component.hatch == False and component.fill:
                    component.edge_color = "None"

                component.fill_plot = ax_hist.fill_between(hx, np.zeros(len(hy)), hy,
                facecolor=component.face_color, edgecolor = component.edge_color,
                alpha=component.fill_alpha, hatch=component.hatch, **component.fill_kwargs)
                legend_handles.append( (component.plot[0], component.fill_plot) )
            else:
                legend_handles.append(component.plot[0])

            legend_labels.append(component.title)

    #total pdf - essentially a copy paste from pdf.get_curve
    total_h = model.createHistogram(observable.GetName(), npoints_curve)
    if total_h.GetNbinsX() != npoints_curve: #model must have intrinsic binning
        total_hy, total_hx_bin_edges = th12uarray(total_h)
        total_hy     = unp.nominal_values(total_hy)
        total_hy    /= total_hy.sum()
        total_hy    *= data_uarray.sum().n
        total_hx     = (total_hx_bin_edges[:-1] + total_hx_bin_edges[1:])/2

        total_fit_plot = ax_hist.hist(total_hx, weights=total_hy,
                                      bins=total_hx_bin_edges, color=fit_color, histtype='step')
        legend_handles.append(total_fit_plot[2][0])

    else:
        h = model.createHistogram(observable.GetName(), npoints_curve)
        total_hy, total_hx = root_numpy.hist2array(h, False, False, True)
        total_hy = npoints_curve * total_hy / np.sum(total_hy)
        total_hx = (total_hx[0][:-1] + total_hx[0][1:])/2.

        #normalise to the sum of events in plot as per the roofit manual page 12
        #Note that the normalization of the PDF, which has an intrinsic normalization to unity by definition, is automatically adjusted to the number of events in the plot.
        total_hy = np.array(total_hy * data_uarray.sum().n )/ nbins
        total_fit_plot = ax_hist.plot(total_hx, total_hy, color=fit_color)
        legend_handles.append(total_fit_plot[0])

    legend_labels.append(legend_fit_name)

    if legend:
        ax_hist.legend(legend_handles, legend_labels, **legend_kwargs)

    #remove the extra space pyplot adds to the side
    ax_hist.set_ylim(bottom=0)
    ax_hist.set_xlim(observable.getMin(), observable.getMax())

    #calculate and plot the pull
    pull = (unp.nominal_values(data_uarray) - np.interp(bin_centers, total_hx, total_hy)) / unp.std_devs(data_uarray)

    fillBetweenBins = []
    for i in range(len(bins)):
        if ((i == 0) or (i==len(bins)-1)):
            fillBetweenBins.append(bins[i])
        else:
            fillBetweenBins.append(bins[i])
            fillBetweenBins.append(bins[i])

    plt.axhline(0,color='gray', ls=':')
    fillBetweenPlotSet = []
    for i in range(len(pull)):
        fillBetweenPlotSet.append(pull[i])
        fillBetweenPlotSet.append(pull[i])

    ax_pull.fill_between(fillBetweenBins, 0, fillBetweenPlotSet, facecolor='lightgray', edgecolor='lightgray', alpha=1)
    ax_pull.errorbar(bin_centers, pull, xerr= bin_widths/2, fmt='.',color='Black', lw=1.0, capthick = 0)#, xerr=binHalfWidths)
    ax_pull.grid()

    pull_low, pull_up = ax_pull.get_ylim()
    if pull_low > -2: pull_low = -2.1
    if pull_up  <  2: pull_up  =  2.1
    ax_pull.set_ylim(pull_low, pull_up)

    # Titles and labels
    ax_hist.set_title(title, **title_kwargs)
    ax_pull.set_xlabel("%s [%s]" % (observable.GetTitle(),observable.getUnit()), **xlabel_kwargs)#, ha='right', x=1.0)
    ylabel_pull = ax_pull.set_ylabel('Pull', *pull_ylabel_kwargs)
    ylabel_hist = ax_hist.set_ylabel('Events / (%.3f %s)' % (min(bin_widths),observable.getUnit()), **ylabel_kwargs)

    #align the y labels - crude but works for now. align_ylabels only works with gridspec
    plt.gcf().canvas.draw()
    new_label_xcoord = np.min([x for x,y in [ylabel_hist.get_position(), ylabel_pull.get_position()]])
    fig_size = fig.get_size_inches()*fig.dpi # size in pixels
    new_label_xcoord /= fig_size[0] #transform to fractional coordinates
    ax_hist.yaxis.set_label_coords(new_label_xcoord,ylabel_hist.get_position()[1])
    ax_pull.yaxis.set_label_coords(new_label_xcoord,ylabel_pull.get_position()[1])

    if filename != None:
        plt.savefig(filename, bbox_inches='tight')
    if show:
        plt.show()

    return
