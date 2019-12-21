// 平均 mu、標準偏差 sigma の正規分布を生成し、その確率分布を求め、図示する。
// 一様乱数の和を用いる方法による
// (sigma > 0)
clear

// 乱数の生成
function x = normrand2(mu,sigma)
    // ***未完***
    n_sum = 100.0
    sum = 0.0
    for k = 1:n_sum
      sum = sum + rand()
    end
    y = sum / n_sum
    x = 2.0 * sqrt(3.0 * n_sum) * sigma * (y - 1.0 / 2.0) + mu
endfunction

// 以下、この関数を用いたテストを行う。
n = 10000; // 生成する乱数の個数
mu = 5.0; // 乱数のパラメタ
sigma = 1.0;

fmin = 0; fmax = 10;  // 度数分布を求める x の範囲
fn = 40; // その範囲の区間数
df = (fmax-fmin)/fn; // 各区間の幅
fx = fmin + (0:fn)*df;  // 区間の境界の値（両端を含む）
frequency = zeros(1,fn); // 度数分布を納める横ベクトル

an = 100; // 累積確率分布を求める際の点数
da = (fmax-fmin)/an; // 各区間の幅
ax = fmin + (0:an)*da;  // 区間の境界の値（両端を含む）
accumulate = zeros(1,an); // 度数分布を納める横ベクトル

// n 個の乱数を生成し、度数分布と累積度数分布を求める。
// ここでは x < fmin または x >= fmax となる x は数えずに捨てることにした。
for i=1:n
  x = normrand2(mu,sigma);
  if x >= fx(1) then  // x < fx(1) = fmin でないときに、以下を行う。
    // 度数分布では、どの区間に入るか判定し、そこの度数を 1 増やす。
    for j=1:fn
      if x < fx(j+1) then
        frequency(j) = frequency(j)+1;
        break;
      end
    end
    // 累積度数分布では、どの区間に入るか判定し、その区間以降の度数を増やす。
    for j=1:an
      if x < ax(j+1) then
        accumulate(j:an) = accumulate(j:an) + ones(1,an-j+1);
        break;
      end
    end
  end
end

// 結果の画面への表示（不要のときにはコメントアウト）
//disp(frequency);
//disp(accumulate);

// の図示
// 横軸は区間の中心（左右の境界の平均値）をとることにする。
pfx = ( fx(1:fn) + fx(2:fn+1) )/2;
pax = ( ax(1:an) + ax(2:an+1) )/2;
clf;

// 度数分布と累積度数分布を確率と累積確率分布に直して図示する。
frequency = frequency/n;
accumlate = accumulate/n;

clf;
subplot(2,1,1);  // 確率分布 - さきほどの度数分布と同じ。
maxfreq = max(frequency);
plot(pfx,frequency,'bo');
plot([pfx; pfx],[frequency; zeros(frequency)],'b-');
a=gca();
a.data_bounds=[fmin -maxfreq*0.05;fmax maxfreq*1.05];
xlabel('x'); ylabel('p(x)');
title('生起確率');

subplot(2,1,2);  // 累積確率分布
//plot(pax,accumlate,'b.');
plot(pax,accumlate,'b-');
plot(pax,(1+erf((pax-mu)/sqrt(2)/sigma))/2,'r-'); // 累積分布の理論値
a=gca();
a.data_bounds=[fmin -0.05;fmax 1.05];
xlabel('x'); ylabel('F(x)');
title('累積確率分布');
