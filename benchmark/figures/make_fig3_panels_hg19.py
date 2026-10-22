# Supplementary Figures: panels A and B of Fig. 3 rendered separately for GRCh37.
import sys, matplotlib
src=open('make_fig3_combined.py',encoding='utf-8').read().split('for build in BUILD_ORDER:')[0]
exec(src)
out=sys.argv[1]
for name,fn,size in (('category',panel_a,(20,12.5)),('manual',panel_b,(20,16))):
    fig=plt.figure(figsize=size, layout="constrained")
    fn(fig,'hg19')
    for ext in ('pdf','png'):
        fig.savefig(f"{out}/fig3_{name}_hg19.{ext}", dpi=150, bbox_inches="tight")
    plt.close(fig)
print('ok')
