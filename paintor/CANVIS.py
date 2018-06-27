# adaptation of PAINTOR CANVIS.py script -- https://github.com/gkichaev/PAINTOR_V3.0

from optparse import OptionParser
import sys
import subprocess
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
from scipy.stats import norm
import math
import svgutils.transform as sg
import warnings
import os
from reportlab.graphics import renderPDF, renderPM
from reportlab.platypus import SimpleDocTemplate, Image, Indenter
from svglib.svglib import svg2rlg
from reportlab.lib.pagesizes import letter, A5, inch

class RotatedImage(Image):
    def wrap(self,availWidth,availHeight):
        h, w = Image.wrap(self,availHeight,availWidth)
        return w, h
    def draw(self):
        self.canv.rotate(45)
        Image.draw(self)

def vararg_callback(option, opt_str, value, parser):
    """Function that allows for a variable number of arguments at the command line"""
    assert value is None
    value = []
    def floatable(str):
        try:
            float(str)
            return True
        except ValueError:
            return False

    for arg in parser.rargs:
        if arg[:2] == "--" and len(arg) > 2:
            break
        if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):
            break
        value.append(arg)
    del parser.rargs[:len(value)]
    setattr(parser.values, option.dest, value)

def Read_Input(locus_file, data_names, ld_files, annotation_file, annotation_cols, threshold = np.inf):
    """Function that reads in all your data files"""
    # read data and subset by threshold
    paintor_data = pd.read_csv(locus_file, delim_whitespace=True)
    paintor_data = paintor_data[paintor_data[data_names[0]] < threshold]
    old_index = paintor_data.index.values
    paintor_data = paintor_data.reset_index(drop=True)

    data_toplot = paintor_data[data_names]
    position = paintor_data['pos']
    posterior = paintor_data['Posterior_Prob']

    lds = []
    if ld_files is not None:
        for ld_file in ld_files:
            ld = pd.read_csv(ld_file, header=None, delim_whitespace=True)
            ld = ld.iloc[old_index,old_index]
            ld = ld.reset_index(drop=True)
            ld.columns = range(ld.shape[1])
            lds.append(ld)
    else:
        lds = None
    if annotation_file is not None:
        annotation_data = pd.read_csv(annotation_file, delim_whitespace=True)
        if annotation_cols is not None:
            annotation_data = annotation_data[annotation_cols]
        else: 
            header = pd.read_csv(annotation_file, delim_whitespace=True, header=None)
            header = header.values.tolist()
            annotation_cols = header[0]
            annotation_data = annotation_data[annotation_cols]
        annotation_data = annotation_data.loc[old_index,:]
        annotation_data = annotation_data.reset_index(drop=True)
        annotation_data = annotation_data.values
    else:
        annotation_data = None

    data_toplot = data_toplot.values
    position = position.values
    posterior = posterior.values


    return [data_toplot, posterior, position, lds, annotation_data, annotation_cols]

def Plot_Position_Value(position, posterior, threshold, greyscale):
    """Function that plots z-scores, posterior probabilites, other features """
    # make the plot greyscale if desired
    if greyscale == "y":
        plot_color = '#BEBEBE'
        set_color = '#000000'
    else:
        plot_color = '#2980b9'
        set_color = '#D91E18'
    # get the credible set
    [credible_loc, credible_prob] = Credible_Set(position, posterior, threshold)
    # start figure
    fig = plt.figure(figsize=(6, 3.25))
    # subplot for posterior plot
    sub1 = fig.add_subplot(1,1,1, facecolor='white')
    plt.xlim(np.amin(position), np.amax(position)+1)
    plt.ylabel('Posterior probabilities', fontsize=10)
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.xlabel('Position', fontsize=10)
    sub1.scatter(position, posterior, color=plot_color, clip_on=False)
    if threshold != 0:
        sub1.scatter(credible_loc, credible_prob, color=set_color, label='Credible Set', clip_on=False)
        title = "Credible Set: " + str(threshold*100) + "%"
        credible_set = mpatches.Patch(color=set_color, label=title)
        legend = plt.legend(handles=[credible_set])
        for label in legend.get_texts():
            label.set_fontsize(10)
    plt.gca().set_ylim(bottom=0)
    posterior_plots = fig
    return posterior_plots

def Credible_Set(position, posterior, threshold):
    """Function that finds the credible set according to a set threshold"""
    total = sum(posterior)
    bounds = threshold*total
    #make into tuples
    tuple_vec = []
    for i in range(0, len(position)):
        tup = (position[i], posterior[i])
        tuple_vec.append(tup)
    #order tuple from largest to smallest
    tuple_vec = sorted(tuple_vec, key=lambda x: x[1], reverse=True)
    credible_set_value = []
    credible_set_loc = []
    total = 0
    for tup in tuple_vec:
        total += tup[1]
        credible_set_loc.append(tup[0])
        credible_set_value.append(tup[1])
        if total > bounds:
            break
    return credible_set_loc, credible_set_value

def Plot_Annotations(annotation_cols, annotation_data, greyscale):
    """Plot the annotations with labels"""
    annotation_tuple = []
    for i in range(0, len(annotation_cols)):
        annotation = annotation_data[:,i]
        colors = []
        if greyscale == "y":
            for a in annotation:
                if a == 1:
                    colors.append('#000000')
                else:
                    colors.append('#FFFFFF')
        else:
            color_array = ['#2980b9']
            for a in annotation:
                if a == 1:
                    colors.append(color_array[0])
                else:
                    colors.append('#FFFFFF')
        fig = plt.figure(figsize=(5, .75))
        ax2 = fig.add_axes([0.05, 0.8, 0.9, 0.15])
        cmap = mpl.colors.ListedColormap(colors)
        cmap.set_over('0.25')
        cmap.set_under('0.75')
        bounds = range(1, len(annotation)+1)
        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
        annotation_plot = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional',
                                                    orientation='horizontal')
        annotation_plot.set_label(annotation_cols[i], fontsize=8)
        annotation_plot.set_ticks([])
        annotation_plot = fig
        annotation_tuple.append(annotation_plot)
    return annotation_tuple

def Plot_Statistic_Value(position, data_toplot, values, greyscale, lds, pval):
    """function that plots pvalues from given zscores"""

    values_tuple = []
    for i in range(0, len(values)):
        fig = plt.figure(figsize=(6, 3.25))
        sub = fig.add_subplot(1,1,1, facecolor='white')
        plt.xlim(np.amin(position), np.amax(position) + 1)
        plt.tick_params(axis='both', which='major', labelsize=10)
        plt.ylabel('-log10(pvalue)', fontsize=10)
        v = data_toplot[:, i]

        if pval:
            pvalue = [-math.log(zv,10) for zv in v]
        else:
            pvalue = Zscore_to_Pvalue(v)

        if lds is not None:
            if i < len(lds): # exists a corresponding LD
                correlation_matrix = lds[i]
                [top_vect, top_SNP] = Find_Top_SNP(v, correlation_matrix, pval)

            else: # no corresponding LD, so use previously calculated one
                # warnings.warn("Warning: no corresponding LD matrix for data_toplot. Plot is made using previous LD matrix.")
                n = len(lds) - 1
                correlation_matrix = lds[n]
                [top_vect, top_SNP] = Find_Top_SNP(v, correlation_matrix, pval)

            if greyscale == 'y':
                sub.scatter(position, pvalue, c=top_vect, cmap='Greys', zorder=1, clip_on=False)
            else:
                sub.scatter(position, pvalue, c=top_vect, cmap='jet', zorder=1, clip_on=False)
            x = position[top_SNP]
            y = pvalue[top_SNP]
            sub.plot(x, y, marker='D', color='black', zorder=2, clip_on=False)
        else:
            if greyscale == "y":
                sub.scatter(position, pvalue, color='#6B6B6B', clip_on=False)
            else:
                sub.scatter(position, pvalue, color='#D64541', clip_on=False)

        plt.gca().set_ylim(bottom=0)
        #add threshold line at 5*10-8
        x = [np.amin(position), np.amax(position) + 1]
        y = [-1*math.log(5*10**-8)/(math.log(10)), -1*math.log(5*10**-8)/(math.log(10))]
        plt.plot(x,y,'gray', linestyle='dashed', clip_on=False)
        label = mpatches.Patch(color='#FFFFFF', label=values[i])
        legend = plt.legend(handles=[label])
        for label in legend.get_texts():
            label.set_fontsize('small')
        value_plot = fig

        bar = plt.figure()

        if lds is not None:
            # add color bar
            min_value = np.amin(top_vect)
            max_value = np.amax(top_vect)
            fig = plt.figure(figsize=(3, 1.0))
            ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
            if greyscale == 'y':
                cmap = mpl.cm.binary
            else:
                cmap = mpl.cm.jet
            norm = mpl.colors.Normalize(vmin=min_value, vmax=max_value)
            mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='horizontal')
            bar = fig

        values_tuple.append((value_plot,bar))

    return values_tuple

def Zscore_to_Pvalue(zscore):
    """Function that converts zscores to pvalues"""
    abs_zscore = np.absolute(zscore)
    pvalue = -1 * (norm.logsf(abs_zscore) / math.log(10))
    return pvalue

def Find_Top_SNP(value_vect, correlation_matrix, pval):
    correlation_matrix = correlation_matrix.values
    # use r^2
    correlation_matrix = np.square(correlation_matrix)
    value_vect = np.absolute(value_vect)
    if pval:
        top_SNP = value_vect.argmin()
    else:
        top_SNP = value_vect.argmax() # returns index
    # get column corresponding to top SNP
    top_vect = correlation_matrix[:][top_SNP]
    return top_vect, top_SNP

def Plot_Heatmap(lds, greyscale, large_ld):
    """Function that plots heatmap of LD matrix"""
    ld_arr = []
    for correlation_matrix in lds:
        n = correlation_matrix.shape
        correlation_matrix = np.square(correlation_matrix)
        if n[0] > 350 and large_ld == 'n':
            warnings.warn('LD matrix is too large and will not be produced. To override, add "-L y"')
            heatmap = None
            bar = None
            ld_arr.append((heatmap, bar))
            break

        if n[0] > 350 and large_ld == 'y':
            warnings.warn('LD matrix is too large but will be produced due to override flag')

        fig = plt.figure(figsize=(3.25, 3.25))
        sns.set(style="white")
        correlation = correlation_matrix.corr()
        mask = np.zeros_like(correlation, dtype=np.bool)
        mask[np.triu_indices_from(mask)] = True
        if greyscale == "y":
            cmap = sns.light_palette("black", as_cmap=True)
        else:
            cmap = None
        sns.heatmap(correlation, mask=mask, cmap=cmap, square=True,
                    linewidths=0, cbar=False, xticklabels=False, yticklabels=False, ax=None)
        heatmap = fig

        matrix = correlation_matrix.values
        min_value = np.amin(matrix)
        max_value = np.amax(matrix)
        fig = plt.figure(figsize=(3, 1.0))
        ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])
        if greyscale == 'y':
            cmap = mpl.cm.binary
        else:
            cmap = mpl.cm.coolwarm
        norm = mpl.colors.Normalize(vmin=min_value, vmax=max_value)
        mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='horizontal')
        bar = fig
        ld_arr.append((heatmap, bar))
    return ld_arr

def Assemble_PDF(data_plots, posterior_plots, heatmaps, annotation_plot, output):
    DPI = 100
    
    # save posterior plot
    posterior_plots.set_size_inches(8, 4) 
    posterior_plots.savefig('posterior_plots.png', format='png', dpi=DPI, transparent=True, bbox_inches='tight')

    # save annotation plots
    for i,plot in enumerate(annotation_plot):
        plot.set_size_inches(6.9, 0.2)
        plot.savefig('annotation_plot'+str(i)+'.png', format='png', dpi=DPI, transparent=True, bbox_inches='tight')

    # save data plot
    for i,plot in enumerate(data_plots):
        plot[0].set_size_inches(8, 4) 
        plot[0].savefig('value_plots'+str(i)+'.png', format='png', dpi=DPI, transparent=True, bbox_inches='tight')

    # save heatmap
    for i,heatmap in enumerate(heatmaps):
        plot = heatmap[0]
        plot.set_size_inches(np.sqrt(32), np.sqrt(32)) 
        plot.savefig('heatmap'+str(i)+'.png', format='png', dpi=DPI, transparent=True, bbox_inches='tight')

    # colorbar
    colorbar_h = heatmap[1]
    colorbar_h.set_size_inches(4, .5)
    colorbar_h.savefig('colorbar_h.png', format='png', dpi=DPI, transparent=True, bbox_inches='tight')

    colorbar = data_plots[0][1]
    colorbar.set_size_inches(4, .5)
    colorbar.savefig('colorbar.png', format='png', dpi=DPI, transparent=True, bbox_inches='tight')

    # start pdf doc
    doc = SimpleDocTemplate(output+".pdf", pagesize=(10*inch,18*inch+0.5*len(annotation_plot)*inch), rightMargin=0, leftMargin=0, topMargin=0, bottomMargin=0)

    # list for figures
    Story = []

    # add posterior plot
    Story.append(Image('posterior_plots.png'))

    # loop through annotations adding an indent to each one
    for i in range(0,len(annotation_plot)):
        Story.append(Indenter(left=.6*inch))
        Story.append(Image('annotation_plot'+str(i)+'.png'))

    Story.append(Indenter(right=.6*inch))

    # add zscore/pvalue plots
    Story.extend([Image('value_plots'+str(i)+'.png') for i in range(0,len(data_plots))])

    # add colorbar
    Story.append(Image('colorbar_h.png'))

    # add ld plot
    for i in range(0,len(annotation_plot)):
        Story.append(Indenter(left=6.9*inch))
        Story.append(RotatedImage('heatmap'+str(i)+'.png'))

    Story.append(Indenter(right=6.9*inch))

    # add other colorbar
    Story.append(Image('colorbar.png'))

    # build image and save
    doc.build(Story)
    

def Assemble_Figure(data_plots, posterior_plots, heatmaps, annotation_plot, output, horizontal):
    """Assemble everything together and return svg and pdf of final figure"""
    DPI = 300
    size_prob_plot = 215
    size_stat_plot = 225
    size_annotation_plot = 30
    size_heatmap = 200

    if heatmaps == None:
        horizontal = 'n'
    elif len(heatmaps)>1:
        horizontal='y'

    if horizontal == 'y':
        size_width = "9in"
        size_height = '9in'
    else:
        size_width = "5in"
        size_height = '11in'

    fig = sg.SVGFigure(size_width, size_height)
    posterior_plots.savefig('value_plots.svg', format='svg', dpi=DPI, transparent=True)
    posterior_plots = sg.fromfile('value_plots.svg')
    plot1 = posterior_plots.getroot()
    if annotation_plot is not None:
        len_ann_plot = (len(annotation_plot))
    else:
        len_ann_plot = 0
    plot1.moveto(0, 0)
    fig.append(plot1)
    # plot heatmap(s)
    len_annotation_plot = size_annotation_plot * (len_ann_plot + 1)

    if heatmaps is not None:
        heatmap_count = 0
        for heatmap in heatmaps:
            plot4 = heatmap[0]
            plot4.savefig('heatmap.svg', format='svg', dpi=DPI, transparent=True)
            plot4 = sg.fromfile('heatmap.svg')
            plot4 = plot4.getroot()

            #transform and add heatmap figure; must be added first for correct layering
            if horizontal=='y':
                y_scale = len_annotation_plot + size_prob_plot +size_heatmap*heatmap_count+ 75
                plot4.moveto(375,y_scale, scale=1.40)
                plot4.rotate(-45, 0, 0)
                fig.append(plot4)
                x_move = 510

            else:
                y_scale = size_stat_plot*len(data_plots) + (size_heatmap+1) * heatmap_count + len_annotation_plot + size_prob_plot + 110
                plot4.moveto(0,y_scale, scale=1.40)
                plot4.rotate(-45, 0, 0)
                fig.append(plot4)
                x_move=110

            heatmap_count = heatmap_count + 1
        colorbar_h = heatmap[1]
        colorbar_h.savefig('colorbar_h.svg', format='svg', dpi=DPI, transparent=True)
        colorbar_h = sg.fromfile('colorbar_h.svg')
        colorbar_h = colorbar_h.getroot()
        colorbar_h.moveto(x_move, y_scale + size_heatmap)
        fig.append(colorbar_h)


    if annotation_plot is not None:
        # transform and add annotations plots
        index = 0
        for plot in annotation_plot:
            plot.savefig('annotation_plot.svg', format='svg', dpi=DPI, transparent=True)
            plot = sg.fromfile('annotation_plot.svg')
            plot3 = plot.getroot()
            y_move = size_prob_plot + size_annotation_plot * (index + 1)
            plot3.moveto(30, y_move, scale=1.05)
            index += 1
            fig.append(plot3)

    #transform and add zscore plots
    index = 1

    for plot in data_plots:
        plot2 = plot[0]
        plot2.savefig('stats_plot.svg', format='svg', dpi=DPI, transparent=True)
        plot2 = sg.fromfile('stats_plot.svg')
        plot2 = plot2.getroot()
        y_move = size_stat_plot * index + len_annotation_plot
        index += 1
        plot2.moveto(0, y_move)
        fig.append(plot2)


    # extract colorbar
    y_move = size_stat_plot * len(data_plots) + len_annotation_plot + size_prob_plot
    plot = data_plots[0]
    colorbar = plot[1]
    colorbar.savefig('colorbar.svg', format='svg', dpi=DPI, transparent=True)
    colorbar = sg.fromfile('colorbar.svg')
    colorbar = colorbar.getroot()
    colorbar.moveto(100, y_move +40)
    fig.append(colorbar)

    #export final figure as a svg and pdf
    svgfile = output+".svg"
    fig.save(svgfile)

    html_file = open(output+".html",'w+')
    html_str = '<img src="'+output+'.svg" >'
    """
    <img src="canvis.svg" >

    """
    html_file.write(html_str)
    html_file.close()
 
def svgToPdfPng(output):
    img = svg2rlg(output+'.svg')
    doc = SimpleDocTemplate(output+'.pdf',
                            pagesize=(6*inch,13.2*inch),
                            rightMargin=0,
                            leftMargin=0,
                            topmargin=0,
                            bottommargin=0)
    docs = []
    docs.append(img)
 
    doc.build(docs)
    # renderPDF.drawToFile(img, output+'.pdf',autoSize=0)
    # renderPM.drawToFile(img, output+'.png', 'PNG')

def main():

    # Parse the command line data
    parser = OptionParser()
    parser.add_option("-l", "--locus_file", dest="locus_file")
    parser.add_option("-v", "--values", dest="values", action='callback', callback=vararg_callback)
    parser.add_option("-a", "--annotation_file", dest="annotation_file")
    parser.add_option("-c", "--annotation_cols", dest="annotation_cols", action='callback', callback=vararg_callback)
    parser.add_option("-r", "--ld_files", dest="ld_files", action='callback', callback=vararg_callback)
    parser.add_option("-t", "--threshold", dest="threshold", default=0)
    parser.add_option("-g", "--greyscale", dest="greyscale", default='n')
    parser.add_option("-o", "--output", dest="output", default='fig_final')
    parser.add_option("-L", "--large_ld", dest="large_ld", default='n')
    parser.add_option("-H", "--horizontal", dest="horizontal", default='n')
    parser.add_option("-p", "--pval", action='store_true')
    parser.add_option("-T", "--pthresh", dest="pthresh", default=1)

    # extract options
    (options, args) = parser.parse_args()
    locus_file = options.locus_file
    values = options.values
    annotation_file = options.annotation_file
    annotation_cols = options.annotation_cols
    ld_files = options.ld_files
    threshold = options.threshold
    threshold = int(threshold)
    if threshold < 0 or threshold > 100:
        warnings.warn('Specified threshold is not valid; threshold is set to 0')
        threshold = 0
    else:
        threshold = (threshold)*.01
    greyscale = options.greyscale
    output = options.output
    large_ld = options.large_ld
    horizontal = options.horizontal
    pval = options.pval
    pthresh = float(options.pthresh)
    if pthresh < 0:
        warnings.warn('Specified pvalue threshold is not valid; threshold is set to inf')
        pthresh = np.inf

    usage = \
    """ Need the following flags specified (*)
        Usage:
        --locus [-l] specify input file with fine-mapping locus (assumed to be ordered by position)
        --values [-v] specific values to be plotted (column name)
        --annotations [-a]  specify annotation file name
        --annotation_cols [-c] annotations to be plotted
        --ld_files [r] specify the ld matrix file names
        --threshold [-t] threshold for credible set [default: 0]
        --greyscale [-g] sets colorscheme to greyscale [default: n]
        --output [-o] desired name of output file
        --large_ld [-L] overrides to produce large LD despite large size [default: n]
        --horizontal [-h] overrides to produce large LD despite large size [default: n]
        --pval [-p] are the values for plotting pvalues? [default: y]
        --pthresh [-T] threshold over pvalues for plotting [default: inf]
        """

    #check if required flags are presnt
    if(locus_file == None or values == None):
        sys.exit(usage)

    [data_toplot, posterior, position, lds, annotation_data, annotation_cols] = Read_Input(locus_file, values, ld_files, annotation_file, annotation_cols, pthresh)

    # check that we still have some data
    if np.shape(data_toplot)[0] < 1 or np.shape(posterior)[0] < 1 or np.shape(position)[0] < 1 :
        subprocess.call(['touch', output+'.pdf'])
        subprocess.call(['touch', output+'.html'])
        subprocess.call(['touch', output+'.svg'])
        warnings.warn('No variants left to plot, trying removing the pvalue threshold or reviewing the PAINTOR ouput. Returning empty output files.')
        sys.exit(0)

    # plot posterior probability
    posterior_plots = Plot_Position_Value(position, posterior, threshold, greyscale)

    # plot annotations
    if annotation_data is not None:
        annotation_plot = Plot_Annotations(annotation_cols, annotation_data, greyscale)
    else:
        annotation_plot = None

    # plot values vs position, possibly with LD coloring
    data_plots = Plot_Statistic_Value(position, data_toplot, values, greyscale, lds, pval)

    # plot LD
    if lds is not None:
        heatmaps = Plot_Heatmap(lds, greyscale, large_ld)
    else:
        heatmaps = None

    # assemble the whole thing and save
    Assemble_Figure(data_plots, posterior_plots, heatmaps, annotation_plot, output, horizontal)

    Assemble_PDF(data_plots, posterior_plots, heatmaps, annotation_plot, output)    

if __name__ == "__main__":
    main()