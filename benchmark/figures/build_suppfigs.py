# Build the Supplementary Figures PDF. Requires the Arial TrueType fonts (arial.ttf, arialbd.ttf) in ~/.fonts or FONT_DIR.
import io, os, sys
from pypdf import PdfReader, PdfWriter, Transformation, PageObject
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph, Frame
from reportlab.lib.styles import ParagraphStyle
F=os.environ.get('FONT_DIR', os.path.expanduser('~/.fonts'))+'/'; pdfmetrics.registerFont(TTFont('Arial',F+'arial.ttf')); pdfmetrics.registerFont(TTFont('Arial-Bold',F+'arialbd.ttf'))
W,H=595.276,841.89; OUT=sys.argv[2]; SRC=sys.argv[1]
cap=ParagraphStyle('c',fontName='Arial',fontSize=10,leading=13.5,alignment=4)
FIGS=[
 ('criteria_barplot_hg19.pdf','Supplementary Figure S1. Criterion-level performance of ACMG/AMP evidence assignment on the GRCh37 eRepo benchmark (n = 11,563). Sensitivity, specificity, precision, and F1 score (top to bottom) are shown for each of the 18 ACMG/AMP criteria assigned by at least one tool. Performance was compared among InterVar, BIAS-2015, BIAS-2015 with PP5/BP6 disabled, and three HUVar configurations (ESM1b + OpenSpliceAI, AlphaMissense [AM] + SpliceAI, and REVEL/CADD/BayesDel + SpliceAI) using automated evidence only. Criterion assignments were evaluated in a one-versus-rest manner against the corresponding eRepo evidence code. "na" marks criteria that a tool does not evaluate automatically.'),
 ('fig3_category_hg19.pdf','Supplementary Figure S2. Category-level classification performance using automated evidence alone on the GRCh37 eRepo benchmark (n = 11,563). Sensitivity, specificity, precision, and F1 score for pathogenic/likely pathogenic (P/LP), variant of uncertain significance (VUS), and benign/likely benign (B/LB) classifications (top to bottom) are shown, with values indicated above the bars. Performance was compared among InterVar, BIAS-2015, BIAS-2015 without PP5/BP6, and three HUVar configurations (ESM1b + OpenSpliceAI, AlphaMissense [AM] + SpliceAI, and REVEL/CADD/BayesDel + SpliceAI) under both the 2015 ACMG/AMP combining rules (2015) and the Bayesian point-based framework (Bayes).'),
 ('code_ratio_barplot_hg19.pdf','Supplementary Figure S3. Frequency of ACMG/AMP criterion assignment relative to the eRepo reference on the GRCh37 benchmark (n = 11,563). For each criterion, the number of variants assigned the criterion by each tool is shown as the ratio to the number of variants assigned the criterion in the eRepo (ratio 1, black line = same number as the eRepo) for InterVar, BIAS-2015, BIAS-2015 without PP5/BP6, and three HUVar configurations (ESM1b + OpenSpliceAI, AlphaMissense [AM] + SpliceAI, and REVEL/CADD/BayesDel + SpliceAI). Bars above the axis limit are truncated and labelled with their ratio. "na" marks criteria that a tool does not evaluate automatically.'),
 ('fig3_manual_hg19.pdf','Supplementary Figure S4. Effect of incorporating curated evidence into HUVar on the GRCh37 eRepo benchmark (n = 11,563). Sensitivity, specificity, precision, and F1 score are shown for pathogenic/likely pathogenic (P/LP), variant of uncertain significance (VUS), and benign/likely benign (B/LB) classifications across three HUVar configurations (ESM1b + OpenSpliceAI, AlphaMissense [AM] + SpliceAI, and REVEL/CADD/BayesDel + SpliceAI) under both the 2015 ACMG/AMP combining rules (2015) and the Bayesian point-based framework (Bayes). Gray bars indicate performance using automated evidence alone, and blue bars indicate performance after adding the curated evidence of the eRepo supplement (PS2, PS4, PM3, PM6, PP1, PP4, BS3, BS4, BP2, and BP5, and BS2 for genes whose VCEP does not allow gnomAD-based BS2). Values above each pair of bars indicate the performance change.'),
 ('cgp_barplot_hg19.pdf','Supplementary Figure S5. Performance of HUVar in the comprehensive genomic profiling cohort on GRCh37 (n = 1,682 variants). Sensitivity, specificity, precision, and F1 score for pathogenic/likely pathogenic (P/LP), variant of uncertain significance (VUS), and benign/likely benign (B/LB) classification (top to bottom) are shown for three HUVar prediction configurations (ESM1b + OpenSpliceAI, AlphaMissense [AM] + SpliceAI, and REVEL + SpliceAI), for all genes (24 genes) and separately for genes with (12 genes, 1,236 variants) and without (12 genes, 446 variants) VCEP specifications. The reference is the classification of the multidisciplinary expert panel.'),
 ('cgp_sankey3_hg38.pdf','Supplementary Figure S6. Concordance of variant classifications among FoundationOne CDx reports, HUVar, and expert panel curation in the comprehensive genomic profiling cohort (n = 1,682 variants, GRCh38). Sankey diagrams show the F1CDx report, the HUVar Bayesian classification and the multidisciplinary expert panel classification (P/LP, VUS and B/LB) for the three HUVar prediction configurations (rows: ESM1b + OpenSpliceAI, AlphaMissense [AM] + SpliceAI, and REVEL + SpliceAI), for all genes and separately for genes with (12 genes, 1,236 variants) and without (12 genes, 446 variants) VCEP specifications (columns). Panel titles give the proportion of variants for which the HUVar and the F1CDx classifications agree with the expert panel.'),
 ('fig_gene_concordance_revel_manual_hg38_notitle.pdf','Supplementary Figure S7. Agreement of five-tier HUVar classifications with the eRepo by gene on the GRCh38 benchmark (REVEL/CADD/BayesDel + SpliceAI with the curated evidence of the eRepo supplement). For each of the 103 genes with at least 10 eRepo variants, the proportions of variants with exact agreement, a difference in confidence (P vs. LP or B vs. LB), a significant conflict (P/LP or B/LB vs. VUS) and an opposite classification (P/LP vs. B/LB) are shown. Genes are ordered by exact agreement (median 83.9%), and the number of eRepo variants of each gene is given in parentheses.'),
]
def overlay(draw):
    b=io.BytesIO(); c=canvas.Canvas(b,pagesize=(W,H)); draw(c); c.save(); return PdfReader(io.BytesIO(b.getvalue())).pages[0]
w=PdfWriter()
def title(c):
    c.setFont('Arial-Bold',16); 
    for i,t in enumerate(['Supplemental Information for','HUVar: An open framework for automated ACMG/AMP variant classification','incorporating ClinGen Variant Curation Expert Panel specifications']):
        c.drawCentredString(W/2,H-150-24*i,t)
    c.setFont('Arial',13); c.drawCentredString(W/2,H-240,'Hikaru Nakahara et al.')
    for i in range(len(FIGS)): c.drawCentredString(W/2,H-320-22*i,f'Supplementary Figure S{i+1}  Page {i+2}')
    c.setFont('Arial',10); c.drawRightString(W-50,40,'1')
pg=PageObject.create_blank_page(width=W,height=H); pg.merge_page(overlay(title)); w.add_page(pg)
for i,(f,text) in enumerate(FIGS,2):
    fig=PdfReader(os.path.join(SRC,f)).pages[0]; fw,fh=float(fig.mediabox.width),float(fig.mediabox.height)
    p=Paragraph(text,cap); _,ph=p.wrap(W-100,H)
    avail_h=H-60-ph-90; s=min((W-80)/fw,avail_h/fh)
    page=PageObject.create_blank_page(width=W,height=H)
    tx=(W-fw*s)/2; ty=H-40-fh*s
    page.merge_transformed_page(fig,Transformation().scale(s).translate(tx,ty))
    def draw(c,p=p,ty=ty,i=i):
        p.drawOn(c,50,ty-20-ph); c.setFont('Arial',10); c.drawRightString(W-50,40,str(i))
    page.merge_page(overlay(draw)); w.add_page(page)
w.write(OUT); print('written',OUT,len(FIGS)+1,'pages')
