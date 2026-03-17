import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd

# 乱流の影響がある場合のエネルギーカバレッジ確率 理論値
# 乱流の影響がある場合のシミュレーション値

# ============================
# Table I — System Parameters
# ============================
num = 100                   # λ_Lの要素数
H = 5                     # UAVの高さ [m]
D = 0.1                     # 送信機のビーム直径 [m]
Delta_theta = 3.4*10**(-5)  # 送信機のビームの発散角度　[rad]
alpha = 10**(-6)            # 距離減衰係数 [1/m]

p_trans = 600               # レーザーの送信電力 [W]
omega_Achi = 0.004          # ω * A * χ [m^2]
delta_s = 10**(-5)          # エネルギー吸収に割り当てられる受信電力の割合

P_prop = 100                # 固定翼UAVの消費する瞬間的な電力 [W]
Cn2 = 0.5*10**(-14)         # 高度Hにおける屈折率の構造定数　
                            # この定数は10^(-13)~10^(-17)で変化する

# 波数 k=5.92*10^(-3)
k =5.92*10**(6)

# ======================================================
# 収穫電力 P_harv(r)
# ======================================================
def pharv_of_r(r):
    d = np.sqrt(r**2 + H**2)
    P_rec = (omega_Achi * p_trans * np.exp(-alpha * d)) / (D + d * Delta_theta)**2
    return (1 - delta_s) * P_rec 

# ======================================================
# 乱流強度 σ(r)
# ======================================================
def sigma2_r(r):
    return math.sqrt(0.3 * (k**(7/6)) * Cn2 * (r**(11/6))) 

# ======================================================
# Log-normal CDF F_ht(a | σ²(r))
# ======================================================
def F_ht(a, sigma2):
    if a <= 0:
        return 0.0
    mu = -2 * sigma2**2
    s = 2 * sigma2
    return 0.5 * math.erfc(-( (math.log(a) - mu) / (s * math.sqrt(2)) )) 

# ======================================================
# a_r 
# ======================================================
def a_r(r,Pcom):
    return (P_prop+Pcom)/pharv_of_r(r)


# ======================================================
# Energy coverage probability 理論式
# ======================================================
def Penergy_with_turb(lambda_L,P_comm):
    if lambda_L == 0:
        return 0.0
    # 積分区間を等間間隔で区切る
    R_max = 50000
    rs = np.linspace(1e-3, R_max, 2000)
    integrand = []
    for r in rs:
        F = F_ht(a_r(r,P_comm), sigma2_r(r))
        fR = 2 * np.pi * lambda_L * r * np.exp(-np.pi * lambda_L * r**2)
        integrand.append((1 - F) * fR)
    # 台形法による積分
    return np.trapz(integrand,rs)

# ======================================================
# Energy coverage probability 理論式
# ======================================================
def Penergy_with_turb(lambda_L,P_comm):
    if lambda_L == 0:
        return 0.0
    # 積分区間を等間間隔で区切る
    rs = np.linspace(1e-3, 50000, 2000)
    integrand = []
    for r in rs:
        F = F_ht(a_r(r,P_comm), sigma2_r(r))
        fR = 2 * np.pi * lambda_L * r * np.exp(-np.pi * lambda_L * r**2)
        integrand.append((1 - F) * fR)
    # 台形法による積分
    return np.trapz(integrand,rs)

# ======================================================
# Simulation シミュレーション値を計算する
# ======================================================
def Simulation(num_trial, lamL, P_comm):
    if lamL == 0: 
        return 0  # 密度0なら確率は0

    # 成功した回数をカウントする　
    count = 0 

    # 1. 0から1までの一様な乱数を生成 (num_trial個の一括生成)
    u = np.random.uniform(0, 1, num_trial)

    # 2. 最短距離 r を生成
    r_horiz = np.sqrt(-np.log(1 - u) / (np.pi * lamL))

    # --- 修正のポイント：各試行(i)ごとに判定を行う ---
    for i in range(num_trial):
        r = r_horiz[i]
        
        # 3. 乱流係数 ht の生成
        # 対数正規分布の平均(mu)と分散(sigma)は、E[ht]=1となるよう設定するのが一般的
        sig_val = sigma2_r(r)
        mu_val = -0.5 * sig_val**2  # E[ht]=1とするための補正
        # np.random.lognormalは「1つだけ」生成するように修正
        ht = np.random.lognormal(mu_val, sig_val)

        # 4. 受信エネルギーの判定
        # その時のrにおける理想電力に、その時のhtを掛ける
        actual_p = ht * pharv_of_r(r)

        if actual_p >= (P_prop + P_comm):
            count += 1

    return count / num_trial

# ======================================================
# λ_L range
# ======================================================
lambda_vals = np.linspace(0, 1e-6, num)

# P_comm の値ごとに曲線を描く
P_comm_values = [10, 30, 50]
# 理論値を格納する変数
results = {}
# シミュレーション値を格納する変数
result_sim = {}
# シミュレーション回数 100万回
num_trial = 1000000
for P_comm in P_comm_values:
    # 理論値の計算
    results[P_comm] = [Penergy_with_turb(lam,P_comm) for lam in lambda_vals]
    # シミュレーション値の計算
    result_sim[P_comm] = [Simulation(num_trial,lam,P_comm) for lam in lambda_vals]
# ======================================================
# Plot
# ======================================================
plt.figure(figsize=(8,5))
step = 5
for P_comm in P_comm_values:
    # print(P_comm,"[W]","理論値：",results[P_comm],"シミュレーション値:",result_sim[P_comm])

    # 理論値
    plt.plot(lambda_vals, results[P_comm], label=f"Theoretical P_comm = {P_comm} W")
    # シミュレーション値
    plt.plot(lambda_vals[::step],result_sim[P_comm][::step], 
            label=f"Simulated P_comm = {P_comm} W",
            marker="o",markerfacecolor='none',linestyle='none')

plt.xlabel("LBD density λ_L (1/m²)")
plt.ylabel("Energy coverage probability")
plt.grid(True, ls=':')
plt.title("Energy Coverage Probability with Atmospheric Turbulence")
plt.yticks(np.arange(0,1.1,0.1))
plt.ylim(0,1.0)
plt.xlim(0,1e-6)
plt.legend()
plt.show()

# ======================================================
# CSVデータ作成・保存
# ======================================================

# 1. データを格納するための辞書を作成
data_dict = {
    'lambda_L (1/m2)': lambda_vals
}

# 2. P_comm ごとのデータを辞書に追加
for P_comm in P_comm_values:
    # 列名の作成
    theo_col = f'Theoretical_Pcomm_{P_comm}W'
    sim_col = f'Simulated_Pcomm_{P_comm}W'
    
    # リストを辞書に格納
    data_dict[theo_col] = results[P_comm]
    data_dict[sim_col] = result_sim[P_comm]

# 3. pandas DataFrameに変換
df = pd.DataFrame(data_dict)

# 4. CSVファイルとして保存 (Excelで文字化けしないよう utf-8-sig を指定)
file_name = "simulation_results_myver.csv"
df.to_csv(file_name, index=False, encoding='utf-8-sig')

print(f"--- CSV出力完了 ---")
print(f"ファイル '{file_name}' が作成されました。")

# 最初の数行をプレビュー表示
print(df.head())
