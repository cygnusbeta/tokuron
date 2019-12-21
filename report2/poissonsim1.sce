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

function p = expP_theory(t, lambda)
  p = lambda * exp(-lambda * t)
endfunction

function p = logexpP_theory(t, lambda)
  p = -lambda * t + log(lambda)
endfunction

function p = erlangP_theory(k, t, lambda)
  p = (k * lambda) ** k * t ** (k - 1) / factorial(k - 1) .* exp(-k * lambda * t)
endfunction

function p = erlangP_theory_wikipedia(k, t, lambda)
  p = lambda ** k / factorial(k - 1) * t ** (k - 1) .* exp(-lambda * t)
endfunction

function p = logerlangP_theory_wikipedia(k, t, lambda)
  p =  -lambda * t + log(lambda ** k / factorial(k - 1))
endfunction

lambda = 1.0; // (WLOG)単位時間当たり発生頻度 (以降 lambda=1 を仮定)
n=100; // 単位時間の分割数
dt=1.0/n;
p0=lambda*dt; // ~= 各区間での生起確率

kmax=15; tmax=3;  // k, t を求める範囲 (tmax は自然数とする)
poissonP = zeros(tmax,kmax+1);  // Pk(t) を保存する配列
theory = zeros(1,kmax+1); // ある t での理論値を保存する配列

m = 10000; // 試行の数（並列に一挙に全試行を行う。）
events = zeros(m,n*tmax);  // 各試行の、各時間刻みで、事象が生起したか否か。
count = zeros(m,1);  // 事象発生回数のカウンタ
told = zeros(m,1); // 前回事象が発生した時刻
erlangP_len = 6*n*tmax

k_num = 3
k_v = zeros(k_num, 1)

k_v(1) = 1
k_v(2) = 5
k_v(3) = 10
// disp(k_v)
count_k = zeros(k_num, 1)
tgap_k_times = zeros(k_num, 1)
// disp(tgap_k_times)

expP = zeros(n*tmax, 1);
// expP_theory = zeros(n*tmax, 1);
logexpP = zeros(n*tmax, 1);
// logexpP_theory = zeros(n*tmax, 1);
erlangP = zeros(erlangP_len, k_num);
// erlangP_theory_v = zeros(k_num, erlangP_len);

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

      expP(tgap) = expP(tgap) + 1

      // disp('------')
      // disp('i', i)
      // disp('tgap', tgap)
      // disp('tgap_k_times', tgap_k_times)

      // disp(count_k)
      for l = 1:k_num
        count_k(l) = count_k(l) + 1
        tgap_k_times(l) = tgap_k_times(l) + tgap
        if (count_k(l) == k_v(l)) then
          // disp(tgap_k_times(l))
          if (tgap_k_times(l) <= erlangP_len)
            erlangP(tgap_k_times(l), l) = erlangP(tgap_k_times(l), l) + 1
          end
          count_k(l) = 0
          tgap_k_times(l) = 0
        end
      end

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

expP(:) = expP(:)/sum(expP(:))/dt;

t_v_2ntmax = [1:n*tmax]*dt
// disp(t_v_2ntmax)

scf(2); clf;
subplot(3,1,1);
plot(t_v_2ntmax,expP,'ro');
plot(t_v_2ntmax,expP_theory(t_v_2ntmax, lambda),'bo');
title('指数分布');
xlabel('t'); ylabel('P(t)');

logexpP = log(expP)
subplot(3,1,2)
plot(t_v_2ntmax,logexpP,'ro');
// disp('t_v_2ntmax', t_v_2ntmax)
// disp('logexpP_theory(t_v_2ntmax, lambda)', logexpP_theory(t_v_2ntmax, lambda))
plot(t_v_2ntmax,logexpP_theory(t_v_2ntmax, lambda),'bo');
title('指数分布 (片対数グラフ)');
xlabel('t'); ylabel('log(P(t))');

t_v_20ntmax = [1:erlangP_len]*dt

scf(3); clf;
for l = 1:k_num
  subplot(k_num,1,l);
  erlangP_l = zeros(erlangP_len, 1)
  erlangP(:, l) = erlangP(:, l) / sum(erlangP(:, l)) / dt
  erlangP_l = erlangP(:, l)
  // disp(erlangP_l)
  plot(t_v_20ntmax,erlangP_l,'ro');
  // disp('l', l)
  // disp('k_v(l)', k_v(l))
  // disp('lambda', lambda)
  // disp('t_v_20ntmax', t_v_20ntmax)
  plot(t_v_20ntmax, erlangP_theory(k_v(l), t_v_20ntmax, lambda),'bo');
  plot(t_v_20ntmax, erlangP_theory_wikipedia(k_v(l), t_v_20ntmax, lambda),'go');
  plotlabel = 'Erlang 分布 (k = ' + string(k_v(l)) + ')';
  title(plotlabel);
  xlabel('t'); ylabel('P(t)');
  if (l == 1)
    legend(['集計した分布','文献 [2] の定義による理論値','Wikipedia の定義による理論値'], pos='in_upper_right');
  end
end

logerlangP = zeros(erlangP_len, k_num)

scf(4); clf;
for l = 1:k_num
  subplot(k_num,1,l);
  for t = 1:erlangP_len
    logerlangP(t, l) = log(erlangP(t, l) / t ** (k_v(l) - 1))
  end
  logerlangP_l = zeros(erlangP_len, 1)
  logerlangP_l = logerlangP(:, l)
  plot(t_v_20ntmax, logerlangP_l,'ro');
  plot(t_v_20ntmax, logerlangP_theory_wikipedia(k_v(l), t_v_20ntmax, lambda),'go');
  plotlabel = 'Erlang 分布 (f(t)/(t^(k-1)) の片対数グラフ, k = ' + string(k_v(l)) + ')';
  title(plotlabel);
  xlabel('t'); ylabel('log(P(t)/t^(k-1))');
  if (l == 1)
    legend(['集計した分布','Wikipedia の定義による理論値'], pos='in_upper_right');
  end
end

// 結果をコンソールへも表示
// disp(poissonP');

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
