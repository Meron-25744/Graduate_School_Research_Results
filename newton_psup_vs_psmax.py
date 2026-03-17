import numpy as np
import matplotlib.pyplot as plt
import math
import random as rand
# 位置推定した時のレーザーの供給電力に対する給電効率

# 位置推定の試行回数
num=10000
# 基準値eps未満のずれだった割合
Measurement_error_ratio=0
# 試行回数をカウント
trial_count = 0
# 測定誤差
Measurement_error=[[0 for i in range(4)] for j in range(num)]
# 推定座標
Estimated_coordinates=[[0 for i in range(4)] for j in range(num)]
# 真の座標
True_coordinates=[[0 for i in range(4)] for j in range(num)]
# ループ回数
Loop_count=[0 for i in range(num)]
# 実際のドローンの座標と測定した座標の距離の差を格納する
Measurement_diserror = [[0 for i in range(1)] for j in range(num)]

# 部屋の大きさ 単位は[m]
room_X = 5.0
room_Y = 5.0
room_Z = 3.0 

# 最終的な推定値を格納する配列
dX = [0,0,0]


# 各LEDの位置 3っ 単位は[m]
LED_x = [1.0, 4.0, 2.5] 
LED_y = [1.0, 1.0, 4.0]
LED_z = [3.0, 3.0, 3.0]

# 基準値eps 単位は[m]
eps = 0.3

# 今回の関数 二点間の距離の2乗値
def estD(x,y,z,i):
    return (LED_x[i]-x)**2+(LED_y[i]-y)**2+(LED_z[i]-z)**2

# 今回の関数 距離の部分は雑音有と無で変更可
def f(x,y,z,i):
    return math.sqrt(estD(x,y,z,i))-d_true[i]

# 今回使用する関数を偏微分したものの分母
def Rf(x,y,z,i):
    return math.sqrt(estD(x,y,z,i))


# 雑音が入る前の受信電力
pr = [[0 for i in range(3)] for j in range(num)]

# ドローンの位置 単位は[m] ドローンの位置は部屋の範囲内で毎回ランダム
X = rand.uniform(0,5)
Y = rand.uniform(0,5)
Z = rand.uniform(2,3)
# ドローンの座標を変更させる = 初期値(2.5, 2.5, 2.0)とのずれが大きくなり、位置推定の精度が悪くなる

# 推定した各点とドローンの座標との距離の平均
diserror = 0

while trial_count < num:
        
        # 一回のループでの計算回数
        t=0
        
        # 初期値 初期値は部屋の中央付近(2.5, 2.5, 2.0)
        x0 = room_X/2
        y0 = room_Y/2
        z0 = 2.0
        
        # 偏微分項からなるヤコビ行列J
        J = [[0 for i in range(3)] for j in range(3)]
        # 行列Jの逆行列
        J_inv = [[0 for i in range(3)] for j in range(3)]
        # 実際の距離 距離はそれぞれ 単位は[m]
        d_true = [math.sqrt((X-LED_x[0])**2+(Y-LED_y[0])**2+(Z-LED_z[0])**2), 
                  math.sqrt((X-LED_x[1])**2+(Y-LED_y[1])**2+(Z-LED_z[1])**2),
                  math.sqrt((X-LED_x[2])**2+(Y-LED_y[2])**2+(Z-LED_z[2])**2)]
        
        # 各種パラメータ
        Gt = 0  #dBi 送信アンテナ利得
        Gr = 0  #dBi 受信アンテナ利得
        Pt = 28 #W 送信電力 
        r = (3.0*10**8)/(4.5*10**14) #Hz  分母はHz 分子は光の速度m/s  
        d = [0 for i in range(3)]         # 雑音の無い受信電力から計算された距離 真の距離
        SNR = 98 # [dB]  信号対雑音比  送信側

        # 熱雑音電力 = 分散
        N = 28 / 10**(SNR/10) # 単位は[W]
        # 部屋の中央の座標[m]
        Xc=2.5
        Yc=2.5
        Zc=1.5
        # 部屋中央と各LEDとの距離[m]
        d_central = [math.sqrt((Xc-LED_x[0])**2+(Yc-LED_y[0])**2+(Zc-LED_z[0])**2),
                    math.sqrt((Xc-LED_x[1])**2+(Yc-LED_y[1])**2+(Zc-LED_z[1])**2),
                    math.sqrt((Xc-LED_x[2])**2+(Yc-LED_y[2])**2+(Zc-LED_z[2])**2)]
        SNRc_ave =0
        Pc_rW = [0 for i in range(3)]
        SNRc = [0 for i in range(3)]
        # 部屋中央でのSNR 「これを増加させる　=　受信信号の強さ↑　&　雑音電力↓」
        # 部屋中央にドローンがあたった時の無雑音時の平均受信SNR[dB]
        for i in range(3):
            # 受信電力
            Pc_rW[i]= Pt / ( (4*math.pi*d_central[i])/r )**2
            # 受信電力対雑音比
            SNRc[i] = 20 * math.log10(math.sqrt(Pc_rW[i]) / N)
            # SNRの平均
            SNRc_ave += SNRc[i]/3
        
        
        
        # 真の距離から受信電力を計算し、格納するための変数
        Pr_true = [0 for i in range(3)]
        for i in range(3):
            # Pr_true[i] = Pt + Gt + Gr - n*math.log10((4*math.pi*d_true[i])/r)
            Pr_true[i] = Pt / ((4*math.pi*d_true[i])/r)**2
            # 雑音の無い受信電力から計算された距離　確認用
            # d[i] = (r/(4*math.pi))*10**((Pt+Gr+Gt-Pr_true[i])/n)
            # pr[trial_count][i] = 10**(Pr_true[i]/10)
        

        # 確率密度関数　平均値0 雑音はガウス分布に従う np.random.normal(平均,標準偏差,個数)
        Thermal_Noise = np.random.normal(0,N,3)
        
    
        # 受信電力に熱雑音電力を足し合わせる
        Pr_false = [0 for i in range(3)]
        for i in range(3):
            Pr_false[i] = math.sqrt(Pr_true[i]) + Thermal_Noise[i]     # W同士で計算

        

        # 熱雑音を考慮した距離計算
        d_false = [0 for i in range(3)]
        for i in range(3):
            # d_false[i] = (r/(4*math.pi))*10**((Pt+Gr+Gt-Pr_false[i])/n)
            d_false[i] = math.sqrt( ( r/(4*math.pi) )**2 * Pt / Pr_false[i]**2)
        
          
        # ここがニュートン法での計算
        while t < 10:
            Xp = [x0,y0,z0]
            F = [f(x0,y0,z0,0),f(x0,y0,z0,1),f(x0,y0,z0,2)]
            
            J[0][0] = -(LED_x[0]-x0)/Rf(x0,y0,z0,0)
            J[0][1] = -(LED_y[0]-y0)/Rf(x0,y0,z0,0)
            J[0][2] = -(LED_z[0]-z0)/Rf(x0,y0,z0,0)
    
            J[1][0] = -(LED_x[1]-x0)/Rf(x0,y0,z0,1)
            J[1][1] = -(LED_y[1]-y0)/Rf(x0,y0,z0,1)
            J[1][2] = -(LED_z[1]-z0)/Rf(x0,y0,z0,1)
    
            J[2][0] = -(LED_x[2]-x0)/Rf(x0,y0,z0,2)
            J[2][1] = -(LED_y[2]-y0)/Rf(x0,y0,z0,2)
            J[2][2] = -(LED_z[2]-z0)/Rf(x0,y0,z0,2)
            
            # 行列式が0になった時はそもそも計算しない
            if (np.linalg.det(J)==0):
                break
            else:
                J_inv = np.linalg.inv(J)
        
            dX = Xp - np.dot(J_inv,F)
            
            if (np.abs(X-dX[0]) < eps) and (np.abs(Y-dX[1]) < eps) and (np.abs(Z-dX[2]) < eps):
                break
            x0 = dX[0]
            y0 = dX[1]
            z0 = dX[2]
            t+=1
        
        # 結果をまとめて格納
        # ドローンの座標と推定した座標の距離の差を求め、それの平均を求める
        Measurement_diserror[trial_count] = math.sqrt(np.abs(X-dX[0])**2+np.abs(Y-dX[1])**2+np.abs(Z-dX[2])**2)
        
        diserror = Measurement_diserror[trial_count] + diserror
        
        Measurement_error[trial_count] = [np.abs(X-dX[0]),np.abs(Y-dX[1]),np.abs(Z-dX[2]),trial_count]
        Estimated_coordinates[trial_count] = [dX[0],dX[1],dX[2],trial_count]
        True_coordinates[trial_count] = [X,Y,Z,trial_count]
        Loop_count[trial_count] = t
        trial_count+=1



for i in range(num):
    for j in range(3):
        if (Measurement_error[i][j] <= eps) and (Measurement_error[i][j+1] <= eps) and Measurement_error[i][j+2] <= eps:
            # print("受信電力 = ",pr[i],"[W]")
            # print("測定誤差 = ",Measurement_error[i],"[m]")
            # print("推定座標 = ",Estimated_coordinates[i],"[m]")
            # print("真の座標 = ",True_coordinates[i],"[m]")
            # print("ループ回数 = ",Loop_count[i],"[回]")
            # print("各受信機のSNR = ",LED1_SNR[i],LED2_SNR[i],LED3_SNR[i],"[dB]\n")
            Measurement_error_ratio+=1
            break
        else:
            break
        
        
for i in range(num):
    for j in range(3):
        if (Measurement_error[i][j] <= eps) and (Measurement_error[i][j+1] <= eps) and Measurement_error[i][j+2] <= eps:
            break
        else:
            # print("受信電力 = ",pr[i],"[mW]")
            # print("測定誤差 = ",Measurement_error[i],"[m]")
            # print("推定座標 = ",Estimated_coordinates[i],"[m]")
            # print("真の座標 = ",True_coordinates[i],"[m]")
            # print("ループ回数 = ",Loop_count[i],"[回]")
            # print("各受信機のSNR = ",LED1_SNR[i],LED2_SNR[i],LED3_SNR[i],"[dB]\n")
            break

# 結果出力
print("部屋中央での受信機の平均SNR[dB] =",SNRc_ave)
print("送信SNR[dB] =",SNR)
print("\n誤差範囲内に収まった割合 :",Measurement_error_ratio*100/num,"%")
print("基準値 :",eps,"[m]")
print('\n')


#==========================================================================


import matplotlib.pyplot as plt
import numpy as np
import math
# 供給電力に対する最大給電効率

# カーブフィッティング係数
a1 = 0.445
b1 = -0.75
a2 = 0.546
b2 = -0.213

# 各パラメータ
Sigma_1 = 3.92          # 定数
Kappa = 10              #[km] 視界 綺麗な空気なら10
Lambda = 810            #[nm] レーザー波長
Kai = 550               # カイ[nm]
Rho = 1.3               # 光の減衰係数 綺麗な空気なら1.3
Ps = np.arange(1,101,1)                 # 供給電力[w]
# 分散を求める
a=0
for i in range(num):
    a = a + (Measurement_diserror[i]-(diserror/num))**2
# パラメータ設定
radius = 2.5            # 円の半径（cm）
diameter = 2*radius     # 円の直径 (cm)
num_points = 10000      # 点の数
mu = 0                  # 平均

sigma = math.sqrt((a/num)+1)              # 標準偏差 (cm^2)
Parcent_InsideCircle = 0    # 円の中に収まった点の割合
inside_count = 0        # 円の中にある点の数を数える
outside_count = 0       # 円の外にある点の数を数える
e = np.e
ps_num = 100                 # 供給電力の要素数

d = 5.0    # distance [m] 1~10[m] 表示は[m]で
d_2 = d*0.001           # distance [km] 計算は[km]で


trans_ef = [0 for i in range(ps_num)]      # レーザーの伝送効率
distance = [0 for i in range(num_points)]   # ランダムに生成された点の距離

x_inside = [0 for i in range(num_points)]   # 円の中に入った点のx座標
y_inside = [0 for i in range(num_points)]   # 円の中に入った点のy座標

x_outside = [0 for i in range(num_points)]  # 円の外にある点のx座標
y_outside = [0 for i in range(num_points)]  # 円の外にある点のy座標


# ガウス分布から点を生成
# x座標,y座標の長さは原点(0,0)からの距離とする 単位は[cm]
x = np.random.normal(mu, sigma, num_points)
y = np.random.normal(mu, sigma, num_points)


for i in range(num_points):
    distance[i] = math.sqrt(x[i]**2+y[i]**2)
    if distance[i] <= radius:
        x_inside[i] = x[i]
        y_inside[i] = y[i]
        inside_count +=1
    elif distance[i] > radius:
        x_outside[i] = x[i]
        y_outside[i] = y[i]
        outside_count +=1


# x,y_insideとx,y_outside の要素数が過剰なので0の要素を削除
# これで0の要素がグラフに反映されない        
x_inside = [i for i in x_inside if i!=0]
y_inside = [i for i in y_inside if i!=0]
x_outside = [i for i in x_outside if i!=0]
y_outside = [i for i in y_outside if i!=0]

# 点の数をカウントする 
print("Inside Count [個]:",inside_count,"\nOutside Count [個]:",outside_count)
# 円の中の点の個数の割合を求める
Parcent_InsideCircle = 100*(inside_count/num_points)
print("Parcent Inside Circle [%]:",Parcent_InsideCircle)

Eta_max = [0 for i in range(ps_num)]
S_lazer = 3/4    # レーザーの照射面積の割合
Alpha = (Sigma_1/Kappa)*pow((Kai/Lambda),Rho)

# 照射距離に対する最大給電効率,照射距離[m]
for i in range(ps_num):
    trans_ef = (1/e)**(Alpha*d_2)
    Eta_max[i] = S_lazer*(Parcent_InsideCircle/100) * (trans_ef*a2+b2/Ps[i])
    
# 各パラメータを表示
print('半径:[cm]',radius)
print('標準偏差:[cm]',sigma)    
# グラフ描画
plt.title('Power Supply vs Maximum power supply efficiency') 
plt.xlabel('Power Supply [W]')
plt.ylabel('Maximum power supply efficiency')   
plt.plot(Ps,Eta_max)
plt.grid(True)
plt.show()
