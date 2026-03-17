import numpy as np
import matplotlib.pyplot as plt
from scipy.special import lambertw
import matplotlib.patches as pt

# 乱流の影響がないときのシミュレーションと理論値
# コードは一部自分用に改変済み

# ------------------------------ パラメータ ------------------------------
num = 100                               # lamL(x軸)の要素数
H = 5                                   # UAVの高さ [m] 100.0
D = 0.1                                 # 送信機のビーム直径 [m]
theta = 3.4e-5                          # 送信機のビームの発散角度 [rad]
alpha = 1e-6                            # 距離減衰係数

P_trans = 600.0                         # レーザーの送信電力[W]
omega_Achi = 0.004                      # ω * A * χ [m^2]
delta_s = 1e-5                          # エネルギー吸収に割り当てられる受信電力の割合

P_prop = 100.0                          # 固定翼UAVの消費する瞬間的な電力[W]
P_comm_values = [20, 40, 60]            # 通信のために消費する電力[W] 元は10,30,50の3パターン
num_trials = 1000000                    # シミュレーション回数 100万回

result_sim = {}                         # シミュレーションの結果を格納する
lambda_vals = np.linspace(0, 1e-6, num) # レーザー送信機の密度[個/m^2]　変化する
results = {}                            # 理論式の結果を格納する

# ------------------------------ 数式 ------------------------------

# 理論式の中の変数を計算する
def K(P_comm):
    return (alpha/(2*theta))*np.sqrt(((1-delta_s)*omega_Achi*P_trans*np.e**(alpha*D/theta))/(P_prop+P_comm))

# 理論式の中の変数を計算する
def R_star(P_comm):
    return (((2/alpha)*lambertw(K(P_comm), 0).real-D/theta)**2-H**2)**(1/2)

# 理論値を計算する関数
def Penergy_with_turb(lam,P_comm):
    return 1-np.e**(-lam*np.pi*R_star(P_comm)**2)

# 【改善2】simulation関数の引数からcountを削除し、関数内でローカル変数として初期化
# 元のコード: def simulation(num_trials,lamL,Rs,count): 
# 問題点: count引数は外部から渡されるが、関数内でインクリメントしても外部変数には反映されない
#         また、毎回の呼び出しでcountが0から始まることを保証する必要がある
def simulation(num_trials, lamL, Rs):
    count = 0  # 【改善2】関数内で毎回0に初期化することで、各シミュレーションが独立することを保証
    for _ in range(num_trials):
        # 半径R*の円の中に点が何個あるのか
        N = np.random.poisson(lamL * np.pi * Rs**2)
        # 点の個数が1個以上なら、Penergyの定義を満たすと言える
        if N >= 1:
            count += 1
    Penergy_sim = count / num_trials
    return Penergy_sim

# ------------------------------ シミュレーション値と理論値 ------------------------------

# 10,30,50[W]の理論値
for P_comm in P_comm_values:
    arr = []
    for lam in lambda_vals:
        arr.append(Penergy_with_turb(lam, P_comm))
    results[P_comm] = np.array(arr)
    
# 10,30,50[W]のシミュレーション値
for P_comm in P_comm_values:
    arr = []
    Rs = R_star(P_comm)
    for lam in lambda_vals:
        # 【改善3】simulation関数呼び出しからcount引数を削除
        arr.append(simulation(num_trials, lam, Rs))
    result_sim[P_comm] = np.array(arr)
    

# 結果出力
step = 5
for P_comm in P_comm_values:
    print(P_comm,"[W]","理論値：",results[P_comm],"シミュレーション値:",result_sim[P_comm])
    
    plt.plot(lambda_vals, results[P_comm], label=f"P_comm={P_comm} W")
    plt.plot(lambda_vals[::step], result_sim[P_comm][::step], marker=".",markerfacecolor='none',linestyle='None')

plt.grid(True)
plt.yticks(np.arange(0,1.1,0.1))
plt.xlim(0,1.0*10**(-6))
plt.ylim(0,1.0)
plt.show()

