#!/usr/bin/env python3
"""Chinese, source-audited teaching manual. Generates only the requested PDF."""
from pathlib import Path
from io import BytesIO
import math
import re
import json
import os
from xml.sax.saxutils import escape

from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import BaseDocTemplate, PageTemplate, Frame, Paragraph, Spacer, PageBreak, Flowable, Table, TableStyle, KeepTogether
from reportlab.graphics import renderPDF
from matplotlib.mathtext import MathTextParser
from matplotlib.font_manager import FontProperties
from reportlab.graphics.shapes import Drawing, String, Rect

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get('CEE_MANUAL_OUTPUT', str(ROOT/'output/CEE_V0Lite_STAR_KFParticle_Principles_Teaching_Manual_v0.7.0.pdf')))
OUT.parent.mkdir(parents=True,exist_ok=True)
FONT = Path(os.environ.get('CEE_MANUAL_FONT_DIR', str(Path(__file__).resolve().parent/'fonts')))
for name, sub in [('CN','400Regular/NotoSansSC_400Regular.ttf'),('CNB','700Bold/NotoSansSC_700Bold.ttf')]:
    pdfmetrics.registerFont(TTFont(name,str(FONT/sub)))
pdfmetrics.registerFont(TTFont('Mono','/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf'))
pdfmetrics.registerFontFamily('CN',normal='CN',bold='CNB',italic='CN',boldItalic='CNB')
INK=colors.HexColor('#172C3B'); BLUE=colors.HexColor('#12619A'); TEAL=colors.HexColor('#087E79')
GOLD=colors.HexColor('#986017'); PURPLE=colors.HexColor('#715493'); GRAY=colors.HexColor('#51626F')
PALE=colors.HexColor('#F2F6F9'); GRID=colors.HexColor('#C4D2DD'); WHITE=colors.white
W,H=A4; LM=48; CW=W-2*LM
S={
 'body':ParagraphStyle('body',fontName='CN',fontSize=10.8,leading=18.3,textColor=INK,spaceAfter=9,wordWrap='CJK'),
 'h1':ParagraphStyle('h1',fontName='CNB',fontSize=20,leading=28,textColor=INK,spaceAfter=17,wordWrap='CJK'),
 'h2':ParagraphStyle('h2',fontName='CNB',fontSize=12.5,leading=20,textColor=BLUE,spaceBefore=6,spaceAfter=7,wordWrap='CJK'),
 'small':ParagraphStyle('small',fontName='CN',fontSize=9,leading=14.7,textColor=GRAY,spaceAfter=7,wordWrap='CJK'),
 'cap':ParagraphStyle('cap',fontName='CN',fontSize=9,leading=14.7,textColor=GRAY,spaceBefore=5,spaceAfter=11,wordWrap='CJK'),
 'code':ParagraphStyle('code',fontName='Mono',fontSize=8.6,leading=13.4,textColor=INK,spaceAfter=7),
 'cell':ParagraphStyle('cell',fontName='CN',fontSize=9.6,leading=15.5,textColor=INK,wordWrap='CJK'),
 'head':ParagraphStyle('head',fontName='CNB',fontSize=9.6,leading=15.5,textColor=INK,wordWrap='CJK'),
}
def P(s,style='body'):
    s=re.sub(r'<sub>([^<]*)</sub>',lambda m:'<sub>'+m.group(1).replace('ₓ','x').replace('ᵧ','y').replace('p_z','pz')+'</sub>',s)
    for a,b in [('ₓ','<sub>x</sub>'),('ᵧ','<sub>y</sub>'),('ₐ','<sub>a</sub>'),('₀','<sub>0</sub>'),('⁻','<super>−</super>'),('ᵀ','<super>T</super>')]:s=s.replace(a,b)
    s=s.replace('–','-')
    return Paragraph(s,S[style])
def sub(s): return P(s,'h2')
def code(s): return P(escape(s).replace('\n','<br/>').replace('  ','&#160;&#160;'),'code')
def note(title,s):
    t=Table([[P(title,'head')],[P(s,'cell')]],colWidths=[CW-20])
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,-1),PALE),('BOX',(0,0),(-1,-1),.6,GRID),('LEFTPADDING',(0,0),(-1,-1),10),('RIGHTPADDING',(0,0),(-1,-1),10),('TOPPADDING',(0,0),(-1,0),8),('BOTTOMPADDING',(0,-1),(-1,-1),10)]))
    return t
def table(rows,widths=None):
    data=[[P(str(x),'head' if i==0 else 'cell') for x in row] for i,row in enumerate(rows)]
    t=Table(data,colWidths=widths or [CW/len(rows[0])]*len(rows[0]),hAlign='LEFT',repeatRows=1)
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#E5EEF4')),('ROWBACKGROUNDS',(0,1),(-1,-1),[WHITE,colors.HexColor('#F7FAFC')]),('VALIGN',(0,0),(-1,-1),'TOP'),('GRID',(0,0),(-1,-1),.4,GRID),('LEFTPADDING',(0,0),(-1,-1),8),('RIGHTPADDING',(0,0),(-1,-1),8),('TOPPADDING',(0,0),(-1,-1),7),('BOTTOMPADDING',(0,0),(-1,-1),7)]))
    return t

_math_cache={}
_math_parser=MathTextParser('path')
_math_fonts={}
def mdrawing(s,size=14):
    s=re.sub(r'\\(mathbf|mathsf)\s+([A-Za-z])',r'\\\1{\2}',s)
    s=s.replace(r'\mathop{\rm arg\,min}',r'\mathrm{arg\,min}')
    s=re.sub(r'\\frac([0-9])',r'\\frac{\1}',s)
    key=(s,size)
    if key not in _math_cache:
        parsed=_math_parser.parse('$'+s+'$',dpi=72,prop=FontProperties(family='STIXGeneral',size=size))
        d=Drawing(parsed.width,parsed.height+2)
        for font,fs,num,x,y in parsed.glyphs:
            path=font.fname
            if path not in _math_fonts:
                name='MathFont'+str(len(_math_fonts));pdfmetrics.registerFont(TTFont(name,path));_math_fonts[path]=name
            d.add(String(x,y+parsed.depth+1,chr(num),fontName=_math_fonts[path],fontSize=fs,fillColor=INK))
        for x,y,w,h in parsed.rects:d.add(Rect(x,y+parsed.depth+1,w,h,fillColor=INK,strokeColor=None))
        _math_cache[key]=d
    return _math_cache[key]
def drawmath(c,s,x,y,size=13,anchor='left',maxw=None):
    d=mdrawing(s,size); scale=min(1,maxw/d.width) if maxw else 1
    dx=x-(d.width*scale/2 if anchor=='center' else d.width*scale if anchor=='right' else 0)
    c.saveState(); c.translate(dx,y); c.scale(scale,scale); renderPDF.draw(d,c,0,0); c.restoreState()
class Eq(Flowable):
    def __init__(self,s,size=16):
        Flowable.__init__(self); self.s=s; self.d=mdrawing(s,size); self.sc=min(1,(CW-15)/self.d.width); self.width=CW; self.height=self.d.height*self.sc+18
    def draw(self):
        self.canv.saveState(); self.canv.translate((CW-self.d.width*self.sc)/2,10); self.canv.scale(self.sc,self.sc); renderPDF.draw(self.d,self.canv,0,0); self.canv.restoreState()
def eq(s,size=16): return Eq(s,size)
def txt(c,s,x,y,size=10,col=INK,anchor='left',bold=False):
    s=s.replace('⁻','-').replace('–','-')
    c.setFillColor(col); c.setFont('CNB' if bold else 'CN',size)
    if anchor=='center': c.drawCentredString(x,y,s)
    elif anchor=='right': c.drawRightString(x,y,s)
    else:c.drawString(x,y,s)
def para(c,s,x,y,w,size=10):
    st=ParagraphStyle('node',parent=S['body'],fontSize=size,leading=size*1.55,spaceAfter=0)
    p=Paragraph(s,st); _,h=p.wrap(w,1000); p.drawOn(c,x,y-h); return h
def arrow(c,x1,y1,x2,y2,col=BLUE,wd=1.6):
    c.setStrokeColor(col);c.setFillColor(col);c.setLineWidth(wd);c.line(x1,y1,x2,y2)
    a=math.atan2(y2-y1,x2-x1); z=6; p=c.beginPath();p.moveTo(x2,y2)
    p.lineTo(x2-z*math.cos(a-.43),y2-z*math.sin(a-.43));p.lineTo(x2-z*math.cos(a+.43),y2-z*math.sin(a+.43));p.close();c.drawPath(p,fill=1,stroke=0)
class Diagram(Flowable):
    def __init__(self,kind,h=250): self.kind=kind;self.width=CW;self.height=h
    def draw(self):
        c=self.canv;h=self.height
        if self.kind=='cascade':
            pv=(42,98);vx=(220,140);vl=(359,208)
            c.setStrokeColor(BLUE);c.setLineWidth(2.8);p=c.beginPath();p.moveTo(*pv);p.curveTo(90,98,154,121,*vx);c.drawPath(p)
            arrow(c,181,125,*vx,BLUE,2.8);arrow(c,*vx,*vl,TEAL,2.5)
            arrow(c,*vx,346,53,GOLD,2.0);arrow(c,*vl,449,284,PURPLE,2);arrow(c,*vl,465,159,PURPLE,2)
            for x,y,col in [(*pv,BLUE),(*vx,GOLD),(*vl,TEAL)]:c.setFillColor(col);c.circle(x,y,5.5,stroke=0,fill=1)
            txt(c,'碰撞主顶点 PV',14,69,10,BLUE,bold=True);txt(c,'Ξ 衰变点 VΞ',156,165,10,GOLD,bold=True)
            txt(c,'Λ 衰变点 VΛ',298,229,10,TEAL,bold=True)
            drawmath(c,r'\Xi^-\quad S_\Xi',112,134,14,anchor='center');drawmath(c,r'\Lambda\quad S_\Lambda',275,188,14,anchor='center')
            txt(c,'伴随 π⁻',348,39,10,GOLD);txt(c,'质子 p',429,292,10,PURPLE);txt(c,'Λ 的 π⁻',412,142,10,PURPLE)
            txt(c,'PV → VΞ：第一段飞行',14,17,10,BLUE);txt(c,'VΞ → VΛ：第二段飞行',266,17,10,TEAL)
        elif self.kind=='vertex':
            arrow(c,25,54,465,196,BLUE);arrow(c,50,225,462,66,TEAL)
            c.setStrokeColor(GRID);c.setLineWidth(18);c.setStrokeAlpha(.38);c.line(65,67,405,177);c.line(92,209,402,90);c.setStrokeAlpha(1)
            c.setFillColor(GOLD);c.circle(280,136,5,fill=1,stroke=0)
            txt(c,'共同顶点估计',263,164,10,GOLD,bold=True);txt(c,'质子track及其误差范围',31,26,10,BLUE);txt(c,'π track及其误差范围',50,237,10,TEAL)
            txt(c,'真实三维track未必相交；拟合寻找统计上最相容的位置。',15,2,9,GRAY)
        elif self.kind=='poca':
            arrow(c,33,66,464,199,BLUE,2);c.setFillColor(GOLD);c.circle(226,200,5,fill=1,stroke=0)
            c.setFillColor(TEAL);c.circle(247,132,5,fill=1,stroke=0);c.setStrokeColor(GOLD);c.setDash(4,3);c.line(226,200,247,132);c.setDash()
            txt(c,'目标顶点 a',162,218,11,GOLD);txt(c,'track上的最近点 r(s*)',268,112,10,TEAL);txt(c,'动量方向 p',380,210,10,BLUE)
            txt(c,'残差方向垂直于局部切线',40,13,10,GRAY)
        elif self.kind=='constraint':
            c.setFillColor(colors.HexColor('#DCEAF3'));c.ellipse(54,79,138,139,fill=1,stroke=0)
            c.setFillColor(BLUE);c.circle(96,109,4,fill=1,stroke=0)
            c.setFillColor(TEAL);c.circle(402,179,5,fill=1,stroke=0)
            c.setDash(5,3);arrow(c,402,179,56,45,GRAY,1.5);c.setDash();arrow(c,402,179,99,108,BLUE,2.3)
            txt(c,'PV 测量及误差范围',28,156,10,BLUE);txt(c,'Ξ 的衰变区域',342,204,10,TEAL)
            txt(c,'约束前的反向延长线',211,67,10,GRAY);txt(c,'利用 PV 后的track估计',195,151,10,BLUE)
            txt(c,'示意：不是把衰变点挪到 PV，也不是强行穿过 PV 中心。',9,8,9,GRAY)
        elif self.kind=='workflow':
            nodes=[('p 与 Λ 的 π track','输入位置、动量、协方差与粒子假设'),('拟合 Λ 候选','共同顶点 VΛ；合成 Λ 的动量和能量'),('Λ 与伴随 π 组合成 Ξ','把 Λ 当作可输运的复合子粒子，拟合 VΞ'),('Ξ 的拓扑检验','几何筛选；另建副本施加 PV 约束'),('保存、选择并分析','区分未约束候选与各类约束副本')]
            for i,(a,b) in enumerate(nodes):
                y=h-58-i*70;c.setFillColor(PALE);c.setStrokeColor(GRID);c.roundRect(27,y,CW-54,52,6,fill=1,stroke=1)
                txt(c,a,CW/2,y+31,11,BLUE,'center',True);txt(c,b,CW/2,y+12,9,INK,'center')
                if i<4:arrow(c,CW/2,y,CW/2,y-17,GRAY)
        elif self.kind=='endpoints':
            arrow(c,61,107,440,107,TEAL,2)
            for x,label in [(61,'产生端'),(440,'衰变端')]:c.setFillColor(BLUE);c.circle(x,107,5,fill=1,stroke=0);txt(c,label,x,130,11,BLUE,'center',True)
            c.setFillColor(GOLD);c.circle(250,107,5,fill=1,stroke=0);txt(c,'当前参考点',250,83,11,GOLD,'center',True)
            drawmath(c,r'S',250,148,17,anchor='center');txt(c,'完整的产生 → 衰变区间',250,174,10,INK,'center')
            arrow(c,437,57,252,57,GOLD);txt(c,'fSFromDecay：从衰变端到当前点',CW/2,27,10,GOLD,'center')
        elif self.kind=='matrix_uses':
            c.setFillColor(PALE);c.setStrokeColor(GRID);c.roundRect(112,190,272,63,7,fill=1,stroke=1);txt(c,'状态 q 与协方差 C',CW/2,226,13,BLUE,'center',True);txt(c,'中心值 + 误差 + 相互关联',CW/2,205,10,INK,'center')
            labels=[('输运与约束','改变参考点、吸收顶点信息'),('派生物理量','质量、衰变长度、cτ及误差'),('组合更高一级母粒子','Λ 继续参与 Ξ 的顶点拟合')]
            for i,(a,b) in enumerate(labels):
                x=7+i*166;arrow(c,CW/2,190,x+79,139,GRAY);c.setFillColor(PALE);c.setStrokeColor(GRID);c.roundRect(x,35,158,99,5,fill=1,stroke=1)
                para(c,'<b>'+a+'</b>',x+10,115,138,10.2);para(c,b,x+10,82,138,9.6)

class Covariance(Flowable):
    def __init__(self,n):self.n=n;self.width=CW;self.height=(n+1)*32+41
    def draw(self):
        c=self.canv;n=self.n;labels=['x','y','z','p_x','p_y','p_z','E','S'][:n]
        left=48;cell=(CW-57)/n;top=self.height-27
        txt(c,'列变量',CW/2,self.height-7,9,GRAY,'center');txt(c,'行变量',0,top-12,8,GRAY)
        for j,v in enumerate(labels):drawmath(c,v,left+(j+.5)*cell,top-17,12,anchor='center')
        for i,a in enumerate(labels):
            y=top-32*(i+2)
            drawmath(c,a,30,y+10,12,anchor='center')
            for j,b in enumerate(labels):
                x=left+j*cell
                fill=WHITE if i!=j else [colors.HexColor('#E4EFF8')]*3+[colors.HexColor('#E4F2EE')]*3+[colors.HexColor('#EEE8F5'),colors.HexColor('#F8EDD9')]
                c.setFillColor(fill[i] if isinstance(fill,list) else fill);c.setStrokeColor(GRID);c.setLineWidth(.35);c.rect(x,y,cell,32,fill=1,stroke=1)
                s=r'\sigma_{'+a+r'}^2' if i==j else r'C_{'+a+','+b+'}'
                drawmath(c,s,x+cell/2,y+10,10.5,anchor='center',maxw=cell-6)

class TransportMatrix(Flowable):
    def __init__(self):self.width=CW;self.height=210
    def draw(self):
        c=self.canv; labs=['x','y','z','p_x','p_y','p_z'];left=100;cell=51;top=194
        for j,v in enumerate(labs):drawmath(c,v,left+j*cell+25,top,12,anchor='center')
        for i,v in enumerate(labs):
            y=top-30*(i+1);drawmath(c,v+"'",70,y+3,12,anchor='center')
            for j in range(6):
                c.setFillColor(PALE if j==i+3 else WHITE);c.setStrokeColor(GRID);c.rect(left+j*cell,y-3,cell,30,fill=1,stroke=1)
                drawmath(c,'1' if i==j else r'\Delta s' if j==i+3 else '0',left+j*cell+25,y+4,12,anchor='center')
