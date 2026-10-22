#!/usr/bin/env python3
from pathlib import Path
BASE=Path(__file__).resolve().parent
source=(BASE/'rewrite_cee_manual_v070.py').read_text().split("if __name__=='__main__':build()")[0]
exec(compile(source,str(BASE/'rewrite_cee_manual_v070.py'),'exec'),globals())
from reportlab.platypus.tableofcontents import TableOfContents
OLD=PAGES[:];CHAPTERS=[]
def chapter(title,summary):CHAPTERS.append([title,summary,[]])
def section(title,*items,src=''):CHAPTERS[-1][2].append((title,list(items),src))
def reuse(*indices):
    for i in indices:
        title,items,src,_=OLD[i-1];section(title,*items,src=src)
exec((BASE/'book_user_sections.py').read_text(),globals())
chapter('共同基础：状态、协方差与 Kalman 拟合','本章解释两种动量模式共用的几何与误差语言。固定 12D / 7×12 公式是已有两体实现；通用 n-body 扩展输入维数。')
reuse(*range(36,45),*range(51,56),*range(61,70))
chapter('零场 TOF：按执行步骤理解迭代','外层 TOF 自洽迭代更新动量与衰变时刻；内部顶点拟合重线性化。不要把两层循环混为一谈。')
section('先看完整执行顺序',P('① 准备 PV、track 方向与 TOF 匹配；② 直线几何给初始 SV，seedParentBeta 给初始母粒子飞行时间；③ 当前路程与扣除母粒子飞行时间后的 TOF 给 daughter 动量；④ 共同顶点拟合和 Kalman 更新；⑤ 更新母粒子速度和衰变时刻；⑥ 不收敛回到③；⑦ 合成母粒子、传播误差、按需建立 PV 副本。'),P('迭代是因为动量依赖未知 SV 与衰变时刻，而两者又依赖动量。不是重复加入同一个测量来缩小误差。'))
reuse(*range(56,61),70,71,72,81,108)
chapter('有磁 tracking：参考点、POCA 与 SV','动量来自 track 拟合，没有 TOF 外层动量循环；共同顶点内部迭代仍然需要。')
section('先看完整执行顺序',P('① 检查位置、物理动量、电荷和完整协方差；② 默认转换到 track–PV DCA 参考点；③ 直线几何给 SV 初值；④ 原始输入 track 各自求当前 SV 附近 POCA；⑤ 同时更新 SV，不收敛回到④；⑥ 更新 daughter、合成母粒子；⑦ 传播完整拟合误差并选择；⑧ 按配置生成独立 PV 副本并输出。'),P('④–⑤变化的是参考位置、局部线性化与 SV，使用的始终是原测量。当前是初值附近的局部解，遇到多次交叉时不保证遍历所有几何分支。'))
reuse(31,32,33,*range(45,51),74,75)
chapter('PV 约束、S、弧长与寿命','PV 相容性、PV 状态更新和显著性 cut 分别定义。7D 原母粒子与 8D 副本独立保存。')
reuse(77,78,79,80,82,83,84,85,86,87,88,89,91,92,93,94,95,98)
section('级联中的生产点不是都等于碰撞 PV',P('Omega 的生产约束可使用碰撞 PV，Lambda 的生产点则为 Omega SV。Channels 保存 Lambda 原 7D 状态和进入 parent 后的更新副本，但尚未自动为每一级生成相对于各自生产点的独立 8D 飞行状态。不能把 Lambda 对碰撞 PV 的 S/L 当作 Lambda 在 Omega 衰变链里的真实飞行段。'))
chapter('验证、诊断与适用边界','理想几何测试、带误差 toy 与真实实验数据闭合是不同层次的证据。')
section('通道扩展验证记录',table([['检验','结论'],['原核心 / Maker 回归','已有两体、TOF、有磁与接口测试通过'],['两体 / 三体 / 级联','不同 SV、排列顺序、依赖、共享叶 track 排除通过'],['240 个带误差 Lambda toy','240/240；质量 pull 均值 0.1258、宽度 0.9992；平均 χ² 0.9498'],['ROOT 通用树回读','事件号、daughter 链接、角色、协方差与自定义列通过'],['Beno 示例','真实头文件语法检查通过；没有真实 DST 运行闭合'],['实验完整验证','当前没有真实 CEE/FairRoot 环境和样本，不宣称已完成']],[164,CW-164]),P('toy 为受控均匀 By 与高斯误差的两体测试，不代表背景抑制效率，也不代表三体/级联完整 pull 闭合。新增通道没有实测加速倍数承诺。'))
reuse(34,76,90,97,118,119,120)
chapter('附录：CLI、STAR 对照与矩阵索引','旧命令行与两体输出作为兼容参考保留，查阅前先确定当前使用哪条 API。')
reuse(110,111,112,113,114,115,116,117,*range(99,106),121,122,123,124,125,126,127,128)
section('源码与手册重建',P('CEE 源码 SHA256：'+SNAPSHOT_HASH,'small'),P('BenoDst 对照提交：83e58c2f523dcfc1644068afe46c347280e9db80。仅核对上游并提供适配示例，没有修改上游。'),code('CEE_MANUAL_OUTPUT=/absolute/path/manual.pdf \\\n python3 docs/manual_source/build_user_book_v070.py'))
S['h1'].keepWithNext=True;S['h2'].keepWithNext=True
class BookDoc(BaseDocTemplate):
    def __init__(self,path):
        super().__init__(str(path),pagesize=A4,leftMargin=LM,rightMargin=LM,topMargin=60,bottomMargin=58,title='CEE KFParticle 用户与算法手册 v0.7.0',author='CEE KFParticle')
        self.addPageTemplates(PageTemplate(id='book',frames=[Frame(LM,58,CW,H-118,leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)],onPage=self.decorate));self.locations=[]
    def beforeDocument(self):self.locations=[]
    def decorate(self,c,doc):
        c.saveState();txt(c,'CEE KFParticle · 用户与算法手册',LM,H-31,8,GRAY);txt(c,'v0.7.0',W-LM,H-31,8,GRAY,'right');c.setStrokeColor(GRID);c.line(LM,H-41,W-LM,H-41);c.line(LM,40,W-LM,40);txt(c,str(doc.page),W-LM,24,9,GRAY,'right');c.restoreState()
    def afterFlowable(self,f):
        if hasattr(f,'bookkey'):
            self.canv.bookmarkPage(f.bookkey);self.canv.addOutlineEntry(f.booktitle,f.bookkey,f.booklevel,False);self.notify('TOCEntry',(f.booklevel,f.booktitle,self.page,f.bookkey));self.locations.append(dict(title=f.booktitle,page=self.page,level=f.booklevel))
def heading(title,key,level):
    h=P(title,'h1' if level==0 else 'h2');h.booktitle=title;h.bookkey=key;h.booklevel=level;return h
def build_book():
    story=[Spacer(1,90),P('CEE KFParticle<br/>用户与算法手册','h1'),P('v0.7.0 · 自定义衰变通道与分析接口'),Spacer(1,30),P('从输入、编译和用户分析类开始，逐步走到三体、级联、ROOT 输出与重建原理。'),P('输入格式：CEE track/hit 与 BenoDst 适配<br/>重建模式：有磁 tracking 与零场 TOF<br/>状态与输出：原母粒子、PV 副本、daughter 关系及误差'),Spacer(1,40),note('阅读建议','第 1–6 章：实际使用。第 7 章：数据结构与 cut。第 8–11 章：算法原理。第 12 章：验证边界。第 13 章：兼容接口与参考。'),PageBreak(),P('目录','h1')]
    toc=TableOfContents();toc.levelStyles=[ParagraphStyle('toc0',fontName='CNB',fontSize=11,leading=18,spaceBefore=9,textColor=BLUE),ParagraphStyle('toc1',fontName='CN',fontSize=9,leading=14,leftIndent=15,textColor=INK)];story.extend([toc,PageBreak()])
    for ci,(title,summary,sections) in enumerate(CHAPTERS,1):
        if ci>1:story.append(PageBreak())
        story.extend([heading(f'第 {ci} 章　{title}',f'ch{ci}',0),P(summary)])
        for si,(title,items,src) in enumerate(sections,1):
            story.append(heading(f'{ci}.{si}　{title}',f'ch{ci}s{si}',1));story.extend(items)
            if src:story.append(P('源码依据：'+src,'small'))
            story.append(Spacer(1,10))
    doc=BookDoc(OUT);doc.multiBuild(story);(BASE/'book_page_map.json').write_text(json.dumps(doc.locations,ensure_ascii=False,indent=2));print('BOOK',OUT,'chapters',len(CHAPTERS),'sections',sum(len(x[2]) for x in CHAPTERS),'pages',doc.page)
if __name__=='__main__':build_book()
