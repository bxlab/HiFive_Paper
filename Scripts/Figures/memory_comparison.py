#!/usr/bin/env python

import sys
import os
from glob import glob
from math import floor

import numpy
from pyx import canvas, text, path, graph, color, trafo, unit, attr, deco, style, bitmap


unit.set(defaultunit="cm")
text.set(mode="latex")
text.preamble(r"\usepackage{times}")
text.preamble(r"\usepackage{sansmath}")
text.preamble(r"\sansmath")
text.preamble(r"\renewcommand*\familydefault{\sfdefault}")
painter = graph.axis.painter.regular( labeldist=0.1, labelattrs=[text.size(-3)], titleattrs=[text.size(-3)] )
painter2 = graph.axis.painter.regular( labelattrs=None )

step_colors = {
    0:color.cmyk.MidnightBlue,
    1:color.cmyk.Mahogany,
    2:color.cmyk.CadetBlue,
    3:color.cmyk.Dandelion,
    4:color.cmyk.OliveGreen,
}

step_names = {
    '0':'Loading data',
    '1':'Filtering and preprocessing',
    '2':'Binning',
    '3':'Normalization',
    '4':'Heatmapping'
}

def main():
    width = 16.8
    out_fname = sys.argv[1]
    basedir = "%s/Analysis/Timing" % '/'.join(os.path.dirname(os.path.realpath(__file__)).split('/')[:-2])
    data_fnames = {
        "HiFive-Probability":{'0':"%s/hifive_data" % basedir, '1':"%s/hifive_project" % basedir,
                '3':'%s/hifive_prob' % basedir, '4':"%s/hifive_prob_heatmap" % basedir},
        "HiFive-Binning":{'0':"%s/hifive_data" % basedir, '1':"%s/hifive_project_nodist" % basedir,
               '3':'%s/hifive_bin' % basedir, '4':"%s/hifive_bin_heatmap" % basedir},
        "HiFive-Express":{'0':"%s/hifive_data" % basedir, '1':"%s/hifive_project" % basedir,
               '3':'%s/hifive_exp' % basedir, '4':"%s/hifive_exp_heatmap" % basedir},
        "HiFive-ExpressKR":{'0':"%s/hifive_data" % basedir, '1':"%s/hifive_project_nodist" % basedir,
                 '3':'%s/hifive_expKR' % basedir, '4':"%s/hifive_expKR_heatmap" % basedir},
        "HiFive-ExpressKR w/distance":{'0':"%s/hifive_data" % basedir, '1':"%s/hifive_project" % basedir,
                     '3':'%s/hifive_expKRdist' % basedir, '4':"%s/hifive_expKRdist_heatmap" % basedir},
        "HiCPipe":{'0':"%s/bam2raw" % basedir, '1':"%s/hicpipe_data" % basedir,
                   '2':'%s/hicpipe_binning' % basedir, '3':'%s/hicpipe_norm' % basedir,
                   '4':'%s/hicpipe_heatmap' % basedir},
        "HiCLib":{'0':"%s/hiclib_mapping" % basedir, '1':'%s/hiclib_data' % basedir,
                  '3':"%s/hiclib_norm" % basedir, '4':'%s/hiclib_heatmap' % basedir},
        "HiCNorm":{'0':"%s/bam2raw" % basedir, '1':"%s/hicpipe_data" % basedir,
                   "2":"%s/hicnorm_data" % basedir, '3':'%s/hicnorm_norm' % basedir},
    }
    data = load_data(data_fnames)
    c = canvas.canvas()
    c.insert(plot_key(width * 0.3, 1.5), [trafo.translate(width * 0.75 - 1.0, 1.0)])
    c.insert(plot_memory(data, width, 6.0))
    c.writePDFfile(out_fname)

def plot_memory(data, width, height):
    plot_width = width - 4.0
    plot_height = height - 0.8
    left_width = plot_width * 0.6
    right_width = plot_width * 0.2
    prob = data['HiFive-Probability']['3']
    norm = data['HiCNorm']['3']
    prob_min = prob - right_width * 2.2e4 / left_width * 0.7
    prob_max = prob + right_width * 2.2e4 / left_width * 0.3
    norm_min = norm - right_width * 2.2e4 / left_width * 0.6
    norm_max = norm + right_width * 2.2e4 / left_width * 0.4
    c1 = canvas.canvas()
    g1 = graph.graphxy(width=left_width, height=plot_height,
                      y=graph.axis.nestedbar(painter=graph.axis.painter.bar(nameattrs=None)),
                      x=graph.axis.lin(painter=painter, texter=graph.axis.texter.exponential(mantissaexp=r"{{%s}e%s}", nomantissaexp=r"{e%s}"), min=0, max=2.2e4),
                      x2=graph.axis.lin(parter=None, min=0, max=1),
                      y2=graph.axis.lin(painter=None, min=0, max=1))
    c1.insert(g1, [trafo.translate(0, 0)])
    g2 = graph.graphxy(width=right_width, height=plot_height,
                      y=graph.axis.lin(painter=None, min=0, max=1),
                      x=graph.axis.lin(painter=painter, texter=graph.axis.texter.exponential(mantissaexp=r"{{%s}e%s}", nomantissaexp=r"{e%s}"), min=prob_min, max=prob_max),
                      x2=graph.axis.lin(parter=None, min=0, max=1),
                      y2=graph.axis.lin(painter=None, min=0, max=1))
    c1.insert(g2, [trafo.translate(left_width, 0)])
    g3 = graph.graphxy(width=right_width, height=plot_height,
                      y=graph.axis.lin(painter=None, min=0, max=1),
                      x=graph.axis.lin(painter=painter, texter=graph.axis.texter.exponential(mantissaexp=r"{{%s}e%s}", nomantissaexp=r"{e%s}"), parter=graph.axis.parter.linear(tickdists=[5000]), min=norm_min, max=norm_max),
                      x2=graph.axis.lin(parter=None, min=0, max=1),
                      y2=graph.axis.lin(parter=None, min=0, max=1))
    c1.insert(g3, [trafo.translate(left_width + right_width, 0)])
    split = canvas.canvas()
    split.fill(path.path(path.moveto(-0.15, -0.2), path.lineto(0.05, 0.2), path.lineto(.15, 0.2),
               path.lineto(-0.05, -0.2), path.closepath()), [color.cmyk.White])
    split.stroke(path.line(-0.15, -0.2, 0.05, 0.2))
    split.stroke(path.line(-0.05, -0.2, 0.15, 0.2))
    c1.insert(split, [trafo.translate(left_width, 0)])
    c1.insert(split, [trafo.translate(left_width, plot_height)])
    c1.insert(split, [trafo.translate(left_width + right_width, 0)])
    c1.insert(split, [trafo.translate(left_width + right_width, plot_height)])
    methods = ['HiCLib', 'HiCPipe', 'HiCNorm', 'HiFive-Probability', 'HiFive-Binning', 'HiFive-Express',
               'HiFive-ExpressKR', 'HiFive-ExpressKR w/distance']
    hstep = plot_height / len(methods)
    substep = hstep / 6.0
    scale = left_width / 2.2e4
    for i, meth in enumerate(methods[::-1]):
        for j in range(5):
            if str(j) not in data[meth]:
                continue
            if meth not in ['HiCNorm', 'HiFive-Probability'] or j != 3:
                c1.fill(path.rect(0, hstep * i + (4.5 - j) * substep, data[meth][str(j)] * scale, substep),
                       [step_colors[j]])
            elif meth == 'HiCNorm':
                c1.fill(path.rect(0, hstep * i + (4.5 - j) * substep, left_width + right_width * 1.7, substep),
                       [step_colors[j]])
                c1.insert(split, [trafo.translate(left_width, hstep * i + (5 - j) * substep)])
                c1.insert(split, [trafo.translate(left_width + right_width, hstep * i + (5 - j) * substep)])
            else:
                c1.fill(path.rect(0, hstep * i + (4.5 - j) * substep, left_width + right_width * 0.6, substep),
                       [step_colors[j]])
                c1.insert(split, [trafo.translate(left_width, hstep * i + (5 - j) * substep)])
    c = canvas.canvas()
    c.insert(c1, [trafo.translate(4.0, 0.8)])
    for i, meth in enumerate(methods):
        c.text(3.9, height - plot_height / len(methods) * (i + 0.5), meth,
               [text.halign.right, text.valign.middle, text.size(-2)])
    c.text(4.0 + plot_width / 2, 0, "Maximum RAM usage (resident set size, Mbytes)",
           [text.halign.center, text.valign.bottom, text.size(-2)])
    return c

def plot_key(width, height):
    c = canvas.canvas()
    wstep = width / 2.0
    hstep = height / 3
    for i in range(5):
        x = wstep * (i / 3)
        y = height - (i % 3 + 0.5) * hstep
        c.fill(path.rect(x + 0.1, y - 0.15, 0.3, 0.3), [step_colors[i]])
        temp = step_names[str(i)].split(' ')
        if len(temp) > 2:
            c.text(x + 0.5, y + 0.13, ' '.join(temp[:2]),
                   [text.halign.left, text.valign.middle, text.size(-2)])
            c.text(x + 0.5, y - 0.13, temp[-1],
                   [text.halign.left, text.valign.middle, text.size(-2)])
        else:
            c.text(x + 0.5, y, step_names[str(i)],
                   [text.halign.left, text.valign.middle, text.size(-2)])
    return c

def find_max(data):
    max_time = 0.0
    for name in data:
        max_time = max(max_time, find_total(data[name]))
    return max_time

def find_total(data):
    time = 0.0
    for key in data:
        time += data[key]
    return time

def load_data(methods):
    data = {}
    for name in methods:
        data[name] = {}
        for key in methods[name]:
            data[name][key] = find_median(methods[name][key], 'Maximum')
    return data

def find_median(fprefix, name):
    fnames = glob("%s_[0-9].txt" % fprefix)
    if len(fnames) == 0:
        return 0.0
    times = []
    for fname in fnames:
        for line in open(fname):
            temp = line[1:-1].split()
            if len(temp) == 0:
                continue
            if temp[0] == name and name == 'Elapsed':
                temp1 = temp[-1].split(':')
                time = float(temp1[-1])
                for i in range(len(temp1) - 1)[::-1]:
                    time += 60.0 ** (len(temp1) - i - 1) * int(temp1[i])
                times.append(time)
            elif temp[0] == name:
                times.append(float(temp[-1]) / 1000.0)
    if len(times) > 0:
        times = numpy.array(times)
        return numpy.median(times)
    else:
        return 0.0

if __name__ == "__main__":
    main()
