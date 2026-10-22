#!/usr/bin/env python3
from pathlib import Path
import json, hashlib
BASE=Path(__file__).resolve().parent
exec((BASE/'cee_rewrite_layout.py').read_text(),globals())
_paragraph=P
def P(s,style='body'):
    supers=dict(zip('⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾','0123456789+-=()'))
    subs=dict(zip('₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎ᵢⱼₓᵧₛₜ','0123456789+-=()ijxyst'))
    subs['ₖ']='k'
    s=re.sub('['+''.join(supers)+']+',lambda m:'<super>'+''.join(supers[x] for x in m.group())+'</super>',s)
    s=re.sub('['+''.join(subs)+']+',lambda m:'<sub>'+''.join(subs[x] for x in m.group())+'</sub>',s)
    if style!='code':
        for a,b in [('daughter','子粒子'),('mother','母粒子'),('Lambda','Λ'),('pion','π 介子'),('track','track'),('track','track'),('螺旋线','螺旋线')]:
            s=re.sub(r'(?<![\w.\-])'+a+r'(?![\w.\-])',b,s)
    return _paragraph(s,style)

repo=BASE.parents[1]
snapshot=[]
for module in ['CeeV0Lite','CeeV0Analysis']:
    for folder in ['include','src']:
        snapshot += [x for x in (repo/module/folder).glob('*') if x.suffix in ['.h','.cxx']]
for stem in ['CeeTpcTrack','CeeTpcVertex','CeeT0Hit']:
    snapshot += [repo/(stem+ext) for ext in ['.h','.cxx']]
digest=hashlib.sha256()
for f in sorted(snapshot,key=lambda x:x.relative_to(repo).as_posix()):
    digest.update(f.relative_to(repo).as_posix().encode()+b'\0'+f.read_bytes())
SNAPSHOT_HASH=digest.hexdigest()

PAGES=[]
SECTION='导读'
SECTION_RANGES={}
def page(title,*items,src=''):
    PAGES.append((title,list(items),src,SECTION))

class CeeFigure(Flowable):
    def __init__(self,kind,h):self.kind=kind;self.width=CW;self.height=h
    def draw(self):
        c=self.canv
        if self.kind=='decay':
            pv=(48,110);sv=(230,145)
            arrow(c,*pv,*sv,BLUE,2.5);arrow(c,*sv,463,226,TEAL,2.2);arrow(c,*sv,445,58,GOLD,2.2)
            for x,y in [pv,sv]:c.setFillColor(BLUE);c.circle(x,y,5,fill=1,stroke=0)
            txt(c,'碰撞主顶点 PV',13,81,11,BLUE,bold=True);txt(c,'衰变顶点 SV',187,173,11,BLUE,bold=True)
            drawmath(c,r'\Lambda',140,140,21,anchor='center');txt(c,'质子track',375,216,11,TEAL,bold=True);txt(c,'π 介子track',345,72,11,GOLD,bold=True)
            for x,y,col,label in [(463,226,TEAL,'TOF 到达时刻'),(445,58,GOLD,'TOF 到达时刻')]:
                c.setFillColor(col);c.roundRect(x-8,y-8,16,16,3,fill=1,stroke=0);txt(c,label,x-7,y-30,9,col,'right')
            txt(c,'物理示意：先飞行的 Λ，在 SV 处产生两条子粒子track。',9,13,10,GRAY)
        elif self.kind=='loops':
            rows=[('TOF 自洽迭代','给定 SV 和衰变时刻 → 重建动量 → 重新拟合 SV','更新母粒子速度与衰变时刻，再检查是否一致'),('一次顶点拟合','两条原输入track → 各自寻找最近点 → 同时更新顶点','每一轮重新线性化，使用相同原始测量'),('母粒子误差传播','最终 12 维输入 → 每一列取 ± 扰动 → 共 24 次拟合','从输出变化得到 7×12 Jacobian；不重跑完整 TOF 循环')]
            for i,(a,b,d) in enumerate(rows):
                y=self.height-90-i*101;c.setFillColor(PALE);c.setStrokeColor(GRID);c.roundRect(3,y,CW-6,87,5,fill=1,stroke=1)
                txt(c,a,16,y+64,12,BLUE,bold=True);txt(c,b,16,y+40,9.7,INK);txt(c,d,16,y+18,9.3,GRAY)

class ContentsOverview(Flowable):
    def __init__(self):self.width=CW;self.height=420
    def draw(self):
        rows=[['部分','阅读内容／页码']]
        descriptions={'导读':'本书目标与整体计算路线','事件与输入':'PV、track、TOF 匹配与质量假设','状态与误差':'六维、十二维、七维，协方差与导数','重建算法':'输运、POCA、顶点、TOF 闭环与母粒子误差','STAR 对照':'八维状态、产生顶点约束与级联边界','运行与排错':'命令、全部选项、cut、输出与诊断','源码与速查':'文件职责、符号与版本依据'}
        descriptions['八维产生顶点状态']='中性分支、PV 约束、S 与寿命'
        descriptions['用户分析层']='Maker、PID、cut、输出与变量表'
        descriptions['有磁参考点']='默认 DCA、协方差、验证与宏配置'
        descriptions['输运导数优化']='固定路径导数、矩阵计算及对照测试'
        descriptions['带电母粒子扩展']='超氚、By 螺旋线、带电约束与验证'
        for sec,(lo,hi) in SECTION_RANGES.items():rows.append([sec,descriptions.get(sec,'')+'　'+str(lo)+'-'+str(hi)+' 页'])
        t=table(rows,[101,CW-101]);_,h=t.wrap(CW,1000);t.drawOn(self.canv,0,self.height-h)

class MotherCovariance(Covariance):
    def __init__(self):super().__init__(7)
    def draw(self):
        c=self.canv;n=7;labels=['x','y','z','P_x','P_y','P_z','E'];left=48;cell=(CW-57)/n;top=self.height-27
        txt(c,'列变量',CW/2,self.height-7,9,GRAY,'center');txt(c,'行变量',0,top-12,8,GRAY)
        for j,v in enumerate(labels):drawmath(c,v,left+(j+.5)*cell,top-17,12,anchor='center')
        for i,a in enumerate(labels):
            y=top-32*(i+2);drawmath(c,a,30,y+10,12,anchor='center')
            for j,b in enumerate(labels):
                x=left+j*cell;c.setFillColor(PALE if i==j else WHITE);c.setStrokeColor(GRID);c.setLineWidth(.35);c.rect(x,y,cell,32,fill=1,stroke=1)
                s=r'\sigma_{'+a+r'}^2' if i==j else r'C_{'+a+','+b+'}'
                drawmath(c,s,x+cell/2,y+10,10.5,anchor='center',maxw=cell-6)

class LabeledGrid(Flowable):
    def __init__(self,kind):self.kind=kind;self.width=CW;self.height=345 if kind=='joint12' else 240
    def draw(self):
        c=self.canv;labs=['x_1','y_1','z_1','p_{x1}','p_{y1}','p_{z1}','x_2','y_2','z_2','p_{x2}','p_{y2}','p_{z2}'];left=43;cell=(CW-left-4)/12;top=self.height-28;rowh=24
        for j,v in enumerate(labs):drawmath(c,v,left+(j+.5)*cell,top+4,8.5,anchor='center')
        outlabs=labs if self.kind=='joint12' else [r'x_{\rm SV}',r'y_{\rm SV}',r'z_{\rm SV}','P_x','P_y','P_z','E']
        for i,a in enumerate(outlabs):
            y=top-rowh*(i+1);drawmath(c,a,20,y+7,9,anchor='center')
            for j in range(12):
                x=left+j*cell;isblock=(i//6==j//6)
                c.setFillColor(PALE if isblock or self.kind!='joint12' else colors.HexColor('#FFF4EC'));c.setStrokeColor(GRID);c.setLineWidth(.35);c.rect(x,y,cell,rowh,fill=1,stroke=1)
                if self.kind=='joint12':s=('C^{'+str(i//6+1)+'}_{'+str(i%6)+str(j%6)+'}') if isblock else '0'
                elif self.kind=='jacobian':s=r'J_{'+str(i)+','+str(j)+'}'
                else:
                    s='0'
                    if i<3 and (j==i or j==i+6):s=r'\frac{1}{2}'
                    elif 3<=i<6 and (j==i or j==i+6):s='1'
                    elif i==6 and j%6>=3:s=labs[j]+'/E_'+str(j//6+1)
                drawmath(c,s,x+cell/2,y+7,9,anchor='center',maxw=cell-3)
        caption='两个输入 6×6 块独立；交叉块为零是模型假设。' if self.kind=='joint12' else '每一行是一个母粒子输出，每一列是一项子粒子输入。'
        txt(c,caption,8,7,9,GRAY)

class FixedPathMatrix(Flowable):
    def __init__(self):self.width=CW;self.height=245
    def draw(self):
        c=self.canv;labels=['x','y','z','p_x','p_y','p_z'];u=['u_x','u_y','u_z'];left=46;cw=(CW-50)/6;top=220;rh=33
        for j,l in enumerate(labels):drawmath(c,l,left+(j+.5)*cw,top+8,12,anchor='center')
        for i,l in enumerate(labels):
            y=top-(i+1)*rh;drawmath(c,l+"'",22,y+11,12,anchor='center')
            for j in range(6):
                c.setFillColor(PALE if i<3 and j>=3 else WHITE);c.setStrokeColor(GRID);c.rect(left+j*cw,y,cw,rh,fill=1,stroke=1)
                s='1' if i==j else '0'
                if i<3 and j>=3:
                    s=r'\frac{s}{p}('+('1-'+u[i]+'^2' if i==j-3 else '-'+u[i]+u[j-3])+')'
                drawmath(c,s,left+(j+.5)*cw,y+10,10,anchor='center',maxw=cw-7)

class NewtonGraph(Flowable):
    def __init__(self):self.width=CW;self.height=250
    def draw(self):
        c=self.canv;left=50;bottom=39;ww=CW-74;hh=193;xmin=.7;xmax=1.85;ymin=-2;ymax=4.5
        def X(x):return left+(x-xmin)/(xmax-xmin)*ww
        def Y(y):return bottom+(y-ymin)/(ymax-ymin)*hh
        arrow(c,left,Y(0),left+ww+6,Y(0),GRAY,1);arrow(c,X(.75),bottom,X(.75),bottom+hh,GRAY,1)
        path=c.beginPath()
        for i in range(161):
            x=xmin+(xmax-xmin)*i/160;y=x**3-2
            if i==0:path.moveTo(X(x),Y(y))
            else:path.lineTo(X(x),Y(y))
        c.setStrokeColor(BLUE);c.setLineWidth(2);c.drawPath(path)
        sn=1.6;g=sn**3-2;gp=3*sn**2;nxt=sn-g/gp
        c.setStrokeColor(GOLD);c.setLineWidth(1.5);c.line(X(1.14),Y(g+gp*(1.14-sn)),X(1.79),Y(g+gp*(1.79-sn)))
        c.setDash(3,3);c.line(X(sn),Y(0),X(sn),Y(g));c.setDash()
        for x,y,col in [(sn,g,GOLD),(nxt,0,GOLD),(2**(1/3),0,TEAL)]:c.setFillColor(col);c.circle(X(x),Y(y),3,fill=1,stroke=0)
        txt(c,'g(s)=s³−2',X(.91),Y(2.8),11,BLUE);txt(c,'当地切线',X(1.62),Y(3.5),10,GOLD)
        drawmath(c,r's_n',X(sn),Y(0)-24,12,anchor='center');drawmath(c,r's_{n+1}',X(nxt)+7,Y(0)-24,12,anchor='left')
        txt(c,'真根',X(2**(1/3))-2,Y(0)+11,9,TEAL,'right');txt(c,'s',left+ww+5,Y(0)+8,11,GRAY)
        txt(c,'仅用于解释 Newton 的无量纲函数示例；不是候选数据。',9,5,9,GRAY)

def load_part(name,sec):
    global SECTION
    SECTION=sec;exec((BASE/name).read_text(),globals())

load_part('cee_rewrite_foundation_pages.py','导读')
foundation=PAGES[:];PAGES=[]
PAGES.extend(foundation[:4])
load_part('cee_v070_maker_pages.py','用户分析层')
load_part('cee_rewrite_inputs_pages.py','事件与输入')
load_part('cee_v065_reference_pages.py','有磁参考点')
PAGES.extend([(t,it,s,'状态与误差') for t,it,s,_ in foundation[4:]])
load_part('cee_rewrite_math_pages.py','重建算法')
before_extra=len(PAGES)
load_part('cee_rewrite_math_extra_pages.py','重建算法')
extra=PAGES[before_extra:];PAGES=PAGES[:before_extra]
insert_after={
 '先把track移到适合比较的位置':[extra[0]],
 '有磁场时，track为什么变成螺旋线':[extra[1]],
 'Newton 为什么需要二分法保护':[extra[2],extra[3]],
 '从两条 daughter 得到七维母粒子':[extra[4]],
 'TOF 动量为什么离不开衰变顶点':[extra[5]],
 '24 次扰动怎样组成 7×12 导数矩阵':[extra[6]],
}
expanded=[]
for item in PAGES:
    expanded.append(item);expanded.extend(insert_after.get(item[0],[]))
PAGES=expanded
exec((BASE/'validation_content.py').read_text(),globals())
load_part('cee_v064_transport_pages.py','输运导数优化')
load_part('cee_v062_production_pages.py','八维产生顶点状态')
load_part('cee_v063_charged_pages.py','带电母粒子扩展')
load_part('cee_rewrite_comparison_pages.py','STAR 对照')
load_part('cee_rewrite_operations_pages.py','运行与排错')
load_part('cee_rewrite_appendix_pages.py','源码与速查')
PAGES=[(t.replace('daughter','子粒子').replace('pion','π 介子'),it,s,sec) for t,it,s,sec in PAGES]

class SourceMarker(Flowable):
    def __init__(self,s):self.manual_source=s;self.width=0;self.height=0
    def draw(self):pass

class ManualDoc(BaseDocTemplate):
    def __init__(self,path):
        BaseDocTemplate.__init__(self,str(path),pagesize=A4,leftMargin=LM,rightMargin=LM,topMargin=58,bottomMargin=70,title='CEE KFParticle：从探测器输入到 Λ 候选 | v0.7.0 中文修订版',author='CEE KFParticle 中文手册修订')
        f=Frame(LM,70,CW,H-58-70,leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)
        self.addPageTemplates(PageTemplate(id='main',frames=[f],onPage=self.decorate));self.actual=[];self.section=''
    def decorate(self,c,doc):
        c.saveState();txt(c,'CEE KFParticle  |  Λ 重建手册',LM,H-31,8,GRAY);txt(c,'v0.7.0 · 用户分析层',W-LM,H-31,8,GRAY,'right')
        c.setStrokeColor(GRID);c.setLineWidth(.5);c.line(LM,H-40,W-LM,H-40);c.line(LM,57,W-LM,57)
        txt(c,str(doc.page),W-LM,17,9,GRAY,'right');txt(c,'2026-10-22  ·  输入 / 算法 / 误差 / 运行',LM,17,8,GRAY);c.restoreState()
    def afterFlowable(self,f):
        if getattr(f,'manual_title',None):
            title=f.manual_title;key='p'+str(len(self.actual)+1);self.canv.bookmarkPage(key)
            if f.manual_section!=self.section:
                self.section=f.manual_section;self.canv.addOutlineEntry(self.section,key,0,False)
            self.canv.addOutlineEntry(title.replace('\n','：'),key,1,False);self.actual.append({'title':title,'page':self.page,'section':self.section})
        if getattr(f,'manual_source',None):
            st=ParagraphStyle('source',parent=S['small'],fontSize=7.4,leading=10,spaceAfter=0)
            p=Paragraph(f.manual_source.replace('–','-'),st);_,ph=p.wrap(CW-20,100)
            if ph>31:raise ValueError('Source footer too tall: '+f.manual_source)
            p.drawOn(self.canv,LM,50-ph)

def build():
    for i,(_,_,_,sec) in enumerate(PAGES,1):
        if sec not in SECTION_RANGES:SECTION_RANGES[sec]=[i,i]
        SECTION_RANGES[sec][1]=i
    story=[]
    for i,(title,items,src,sec) in enumerate(PAGES):
        if i:story.append(PageBreak())
        label=title if i<4 else f'{i-3:02d}  '+title
        h=P(label.replace('\n','<br/>'),'h1');h.manual_title=title;h.manual_section=sec
        story.append(h);story.extend(items);story.append(SourceMarker(src))
    doc=ManualDoc(OUT);doc.build(story)
    (BASE/'cee_rewrite_page_map.json').write_text(json.dumps(doc.actual,ensure_ascii=False,indent=2))
    print('PDF:',OUT,'topics:',len(PAGES),'pages:',doc.page)
    overflow=[v for i,v in enumerate(doc.actual,1) if v['page']!=i]
    if overflow:print('PAGE_SPILL:',overflow)

if __name__=='__main__':build()
