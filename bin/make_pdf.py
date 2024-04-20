#!/usr/bin/env python

import argparse
import os
from fpdf import FPDF
from datetime import date


argparser = argparse.ArgumentParser(description = 'Combines all generated figures to single PDF. Assumes that figurs are located in the same directory.')
argparser.add_argument('-o', '--output', metavar = 'file', dest = 'out_filename', required = True, help = 'Output file name.')


class MyFPDF(FPDF):        
    def header(self):
        self.set_font('Arial', 'B', 11)
        self.cell(0, 0, f'Report from {date.today().strftime("%B %d, %Y")}', 0, 0, 'L')
        self.ln()
        
    def footer(self):
        self.set_y(-5.5)
        self.set_font('Arial', 'I', 8)
        self.cell(0, 10, 'Page ' + str(self.page_no()), 0, 0, 'C')


if __name__ == '__main__':
    args = argparser.parse_args()

    pdf = MyFPDF(orientation = 'P', unit = 'in', format = 'A4')

    pdf.set_font('Arial', 'B', 15)

    pdf.add_page()
    pdf.image('average_depth.jpeg', x = 1.5, y = 0.7, w = 5, h = 0, type = '', link = '')

    pdf.add_page()
    pdf.image('coverage.jpeg', x = 0.1, y = 0.7, w = 8, h = 0, type = '', link = '')

    pdf.add_page()
    pdf.image('truncated_table.jpeg', x = 0.1, y = 0.7, w = 8, h = 0, type = '', link = '')

    pdf.add_page()
    pdf.image('contamination_boxplot.jpeg', x = 2, y = 0.7, w = 4, h = 0, type = '', link = '')
    pdf.image('contamination_table.jpeg', x = 0, y = 5.7, w = 7, h = 0, type = '', link = '')

    pdf.add_page()
    pdf.image('contamination_per_sample_table.jpeg', x = 0.1, y = 0.7, w = 8, h = 0, type = '', link = '')

    pdf.add_page()
    pdf.image('1000G_pca.jpeg', x = 0.1, y = 0.7, w = 8, h = 0, type = '', link = '')
    pdf.image('1000G_inferred_ancestry_table.jpeg', x = 0.1, y = 8.7, w = 8, h = 0, type = '', link = '')

    pdf.add_page()
    pdf.image('1000G_pca_most_contaminated.jpeg', x = 0.1, y = 0.7, w = 8, h = 0, type = '', link = '')
    
    if os.path.isfile('HGDP_pca.jpeg') and os.path.isfile('HGDP_inferred_ancestry_table.jpeg'):
        pdf.add_page()
        pdf.image('HGDP_pca.jpeg', x = 0.1, y = 0.7, w = 8, h = 0, type = '', link = '')
        pdf.image('HGDP_inferred_ancestry_table.jpeg', x = 0.1, y = 8.7, w = 8, h = 0, type = '', link = '')

    if os.path.isfile('HGDP_pca_most_contaminated.jpeg'):
        pdf.add_page()
        pdf.image('HGDP_pca_most_contaminated.jpeg', x = 0.1, y = 0.7, w = 8, h = 0, type = '', link = '')

    pdf.add_page()
    pdf.image('XY_depth_vs_reported_sex.jpeg', x = 0.1, y = 0.7, w = 8, h = 0, type = '', link = '')

    pdf.output(args.out_filename, 'F')
