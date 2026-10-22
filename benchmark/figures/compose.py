import sys, io
from pypdf import PdfReader, PdfWriter, Transformation, PageObject
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
import glob
import os
ttf=os.path.join(os.environ.get('FONT_DIR', os.path.expanduser('~/.fonts')),'arial.ttf')
pdfmetrics.registerFont(TTFont('Arial',ttf))
W,H=595.276,841.89
def compose(src,label,out,margin=12,top=12):
    fig=PdfReader(src).pages[0]
    fw,fh=float(fig.mediabox.width),float(fig.mediabox.height)
    s=min((W-2*margin)/fw,(H-top-60)/fh)
    page=PageObject.create_blank_page(width=W,height=H)
    tx=(W-fw*s)/2; ty=H-top-fh*s
    page.merge_transformed_page(fig,Transformation().scale(s).translate(tx,ty))
    buf=io.BytesIO(); c=canvas.Canvas(buf,pagesize=(W,H)); c.setFont('Arial',14); c.drawString(14,H-822,label); c.save()
    page.merge_page(PdfReader(io.BytesIO(buf.getvalue())).pages[0])
    w=PdfWriter(); w.add_page(page); w.write(out)

for src,lab,out in [('out/criteria_barplot_hg38.pdf','Fig. 2','final/Fig2.pdf'),('out/fig3_combined_hg38.pdf','Fig. 3','final/Fig3.pdf'),
                    ('out/code_ratio_barplot_hg38.pdf','Fig. 4','final/Fig4.pdf'),('out/cgp_barplot_hg38.pdf','Fig. 5','final/Fig5.pdf'),('out/cgp_sankey3_hg38.pdf','Supplementary Fig. S[X]','final/SuppFig_CGP_sankey.pdf')]:
    compose(src,lab,out)
