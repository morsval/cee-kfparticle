"""Validate derived flight errors with charge/field signs and PV correlations."""
import argparse,json,subprocess
from pathlib import Path
import numpy as np
from equivalent_7d import fit7

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cee-dir',type=Path,required=True);ap.add_argument('--out',type=Path,default=Path('results'));a=ap.parse_args()
    core=a.cee_dir.resolve()/'CeeV0Lite';exe=core/'obj/test_charged_production'
    if not exe.exists():subprocess.run(['make','obj/test_charged_production'],cwd=core,check=True)
    rows=[]
    for idx,line in enumerate(subprocess.check_output([str(exe),'--reference'],text=True).splitlines()):
        f=json.loads(line);y=np.array(f['input']);C=np.array(f['covariance']).reshape(10,10)
        _,_,_,ok,_,L,dl=fit7(y[None,:7],C[None,:7,:7],y[None,7:],C[None,7:,7:],f['by'],f['charge'],C[None,:7,7:],return_flight=True)
        assert ok[0]
        z=np.array(f['state']);c=np.array(f['state_covariance']).reshape(8,8);p=np.linalg.norm(z[3:6]);g=np.zeros(8);g[3:6]=z[7]*z[3:6]/p;g[7]=p
        targetL=z[7]*p;targetdl=np.sqrt(g@c@g)
        rows.append({'id':idx,'charge':f['charge'],'By':f['by'],'cross_covariance_present':bool(np.any(C[:7,7:])),
          'L_absolute_difference_cm':abs(float(L[0])-targetL),'sigmaL_relative_difference':abs(float(dl[0])-targetdl)/targetdl})
    result={'cases':len(rows),'max_L_difference_cm':max(r['L_absolute_difference_cm'] for r in rows),'max_sigmaL_relative_difference':max(r['sigmaL_relative_difference'] for r in rows),'rows':rows}
    result['passed']=bool(result['max_L_difference_cm']<1e-8 and result['max_sigmaL_relative_difference']<1e-8)
    a.out.mkdir(parents=True,exist_ok=True);(a.out/'correlated_fixture_check.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2));assert result['passed']
if __name__=='__main__':main()
