// 再生過程のシミュレーションを行い再生定理を確認する。
clear
xdel(winsid()); // すべてのグラフィックウィンドウを閉じる

// 選択した確率分布での乱数を n 個発生し、横ベクトルで返す。
//   choice = 確率分布の選択番号
//   mu = その分布での期待値
// とりあえず指数分布の場合（choice = 1）しか作っていない。
function p = myrand(n,choice,mu)
  select choice
  case 1 then
    // 平均 mu の指数分布
    p = grand(1,n,'exp',mu);
  else
    disp('unknown choice');
    p = [-n:-1]; // 適当な戻り値を設定した。
  end
endfunction

n=10; // 発生する乱数の個数
mu=0.5; // 乱数の期待値
m=1000; // 図示する際の横軸範囲の分割数
dt=zeros(1,n); // 生起時間の間隔を保存するベクトル
count=zeros(1,m); // 計数関数 N(t) のこと

// 選択した分布での乱数を発生し、生起間隔の時間とする。
dt = myrand(n,1,mu);
// 生起間隔から、生起時刻に換算する。
t(1)=dt(1);
for i=2:n
  t(i)=t(i-1)+dt(i);
end

// 以下に、生起時刻の列から、時刻 t までの生起回数（累計）の関数 N(t) を求め、
// N(t) および N(t)/t を t に対して図示して、初等再生定理と比較してください。
