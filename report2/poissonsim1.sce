// ポアソン分布と指数分布を素過程から求める。
// （指数分布については、求め始めたばかりで未完のプログラム）
// lambda を（WLOG）1 とするとき、k=0,1,2 と t = 1,2,3 での
// Pk(t) を求め図示する。
clear
xdel(winsid()); // すべてのグラフィックウィンドウを閉じる

// 理論的なポアソン分布（比較のため）
function p = poissondist(k,t,lambda)
  // p = 0; // *** とりあえず恒等的に 0 とおいている。正しい式で置き換えてください。
  p = (lambda * t) ** k / factorial(k) * exp(-lambda * t)
endfunction

lambda = 1.0; // (WLOG)単位時間当たり発生頻度 (以降 lambda=1 を仮定)
n=100; // 単位時間の分割数
dt=1.0/n;
p0=lambda*dt; // ~= 各区間での生起確率

kmax=15; tmax=3;  // k, t を求める範囲 (tmax は自然数とする)
poissonP = zeros(tmax,kmax+1);  // Pk(t) を保存する配列
theory = zeros(1,kmax+1); // ある t での理論値を保存する配列

m = 1000; // 試行の数（並列に一挙に全試行を行う。）
events = zeros(m,n*tmax);  // 各試行の、各時間刻みで、事象が生起したか否か。
count = zeros(m,1);  // 事象発生回数のカウンタ
told = zeros(m,1); // 前回事象が発生した時刻

scf(0); clf;
for i=1:n*tmax
  // 確率 p0 の事象が発生したか否かを決める。
  events(:,i) = ( rand(m,1) < p0*ones(m,1) ); // 各メンバーについて T か F が返る。
  // T となった試行について、必要な処理を行う。
  for j=1:m
    if events(j,i) then
      count(j) = count(j)+1;
      // 以下4行は、ポアソン分布を得る際には不要だが、指数分布を得るために加えた。
      tnew = i;  // 発生した区間の番号
      tgap = tnew-told(j);  // 長さ dt の区間数で数えた発生時間間隔
      // *** ここらで時間間隔の頻度を記憶すると良いでしょう。
      told(j) = tnew;
    end
  end
  // t=1,2,... (i=n,2*n,...) ごとに発生回数を集計する。
  if (modulo(i,n)==0) then
    t = i/n;
    for k=0:kmax
      poissonP(t,k+1) = poissonP(t,k+1) + sum(count==k*ones(m,1));
      // count(j)==k となる場合 T となるので、その和を求めると度数が求まる。
      theory(k+1) = poissondist(k,t,lambda);  // 理論値
    end
    // 確率に換算する。
    poissonP(t,:) = poissonP(t,:)/sum(poissonP(t,:));
    // t=1,2,3 では Pk(t) を図示
    if t<=3 then
      subplot(3,1,t);
      plot([0:kmax],poissonP(t,:),'ro');
      plot([0:kmax; 0:kmax],[poissonP(t,:); zeros(poissonP(t,:))],'r-');
      // 理論的な分布と比較
      plot([0:kmax],theory,'b:');
      xlabel('k'); ylabel('Pk(t)');
      plotlabel = 't = ' + string(t);
      title(plotlabel);
      a=gca();
      a.data_bounds=[-0.5 -0.002;5.5 0.502];
    end
  end
end

// 結果をコンソールへも表示
disp(poissonP');

// 最初の３つの試行での事象の発生について別ウィンドウで図示
scf(1); clf;
for j=1:3
  subplot(3,1,j);
  plot([1:n*tmax]*dt,events(j,:));
  xlabel('t'); ylabel('event');
  plotlabel = '試行 ' + string(j);
  title(plotlabel);
  a=gca();
  a.data_bounds=[0 -0.002;tmax 1.002];
end
