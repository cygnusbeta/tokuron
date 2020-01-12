// Markov 連鎖の定常状態を遷移確率行列 P の固有値・固有ベクトルから求める。
// また、連鎖の差分方程式（漸化式）を反復して、定常状態に収束するか調べる。
clear
xdel(winsid()); // すべてのグラフィックウィンドウを閉じる

n = 14

// 状態の数と行列 P を与える。
p = [ 0 0.4 0 0 0 0.3 0 0 0 0 0.3 0 0 0; ...
      0 0.4 0.4 0 0 0 0 0.2 0 0 0 0 0 0; ...
      0.4 0 0 0.3 0 0 0 0 0.3 0 0 0 0 0; ...
      0 0 0 0.8 0.2 0 0 0 0 0 0 0 0 0; ...
      0 0 0 0 0 0 0 0 0 0 0 1 0 0; ...
      0 0 0 0 0 0 1 0 0 0 0 0 0 0; ...
      0 0 0 0 0 0.5 0.5 0 0 0 0 0 0 0; ...
      0 0 0 0 0 0 0 0 1 0 0 0 0 0; ...
      0 0 0 0 0 0 0 0 0 1 0 0 0 0; ...
      0 0 0 0 0 0 0 1 0 0 0 0 0 0; ...
      0 0 0 0 0 0 0 0 0 0 0 0.2 0 0.8; ...
      0 0 0 0 0 0 0 0 0 0 0 0 1 0; ...
      0 0 0 0 0 0 0 0 0 0 0.5 0 0.5 0; ...
      0 0 0 0 0 0 0 0 0 0 0.15 0 0.25 0.6]; // テキスト p.43 の例

p2 = (1 / n) * ones(n, n)
disp('p2')
disp(p2)
d = 0.85
p = (1 - d) * p + d * p2
// n = size(p,1);  // 正方行列 p を与えたと仮定している。
disp('n')
disp(n)
disp('sum(p(i,:)')
for i=1:n
  disp(sum(p(i,:)))
end

// p の左固有ベクトルが欲しい。spec 関数は右固有ベクトルを与えるので、
// p の転置行列の左固有ベクトルを求める。それを転置して、目的である
// p の右固有ベクトルが求まる。
[evec,eval]=spec(p.');
// 対角行列 eval の対角成分が固有値なので取り出す。
lambda = diag(eval);
disp('固有値は');
disp('lambda')
disp(lambda);
// 対応する固有ベクトルが各列に入った行列 evec を転置する。
evec = evec.';
disp('evec')
disp(evec)

// 定常状態の固有ベクトルは、すべての成分が同符号であり、厳密に 1 の
// 固有値に対応していなければならない。それを探す。
for i=1:n
  if abs(sum(evec(i,:)))==sum(abs(evec(i,:))) then // 同符号の成分のとき等号が成立
    if (abs(lambda(i)-1) < 1e-15) then
      disp('定常状態の固有ベクトルが見つかりました。');
      disp([i,lambda(i)]);
      // *** 見つけた evec(i,:) は 2-norm で正規化されています。つまり成分の2乗の和が
      // *** 1 となっています。各成分を定常状態の確率と見るには、1-norm で正規化する
      // *** 必要があります。次の行でそうしてください。
      steadystate = evec(i,:) / sum(evec(i,:));
      disp(steadystate);
      disp('sum(steadystate(i, :))')
      disp(sum(steadystate(1, :)))
      break;
    end
  end
end
if sum(abs(steadystate))==0 then
  halt('定常状態の固有ベクトルが見つかりませんでした。');
end

// 初期条件（時刻 0 での各状態の確率）を適当に与え、行列 p による写像を
// 反復したときの変化を求め、3成分（各状態の確率）がどう変わるか、反復
// 回数を横軸にとって図示してください。
x = [ 0 0 1 0 0 0 0 0 0 0 0 0 0 0 ]; // 初期条件（別のも試してみてください。）
m = 20; // 反復回数

p1 = zeros(1, m + 1)
p2 = zeros(1, m + 1)
p3 = zeros(1, m + 1)
p4 = zeros(1, m + 1)
p5 = zeros(1, m + 1)
p6 = zeros(1, m + 1)
p7 = zeros(1, m + 1)
p8 = zeros(1, m + 1)
p9 = zeros(1, m + 1)
p10 = zeros(1, m + 1)
p11 = zeros(1, m + 1)
p12 = zeros(1, m + 1)
p13 = zeros(1, m + 1)
p14 = zeros(1, m + 1)

p1(1, 1) = x(1, 1)
p2(1, 1) = x(1, 2)
p3(1, 1) = x(1, 3)
p4(1, 1) = x(1, 4)
p5(1, 1) = x(1, 5)
p6(1, 1) = x(1, 6)
p7(1, 1) = x(1, 7)
p8(1, 1) = x(1, 8)
p9(1, 1) = x(1, 9)
p10(1, 1) = x(1, 10)
p11(1, 1) = x(1, 11)
p12(1, 1) = x(1, 12)
p13(1, 1) = x(1, 13)
p14(1, 1) = x(1, 14)

for i=2:m+1
  x = x * p
  p1(1, i) = x(1, 1)
  p2(1, i) = x(1, 2)
  p3(1, i) = x(1, 3)
  p4(1, i) = x(1, 4)
  p5(1, i) = x(1, 5)
  p6(1, i) = x(1, 6)
  p7(1, i) = x(1, 7)
  p8(1, i) = x(1, 8)
  p9(1, i) = x(1, 9)
  p10(1, i) = x(1, 10)
  p11(1, i) = x(1, 11)
  p12(1, i) = x(1, 12)
  p13(1, i) = x(1, 13)
  p14(1, i) = x(1, 14)
end

disp('p1')
disp(p1)
disp('p2')
disp(p2)
disp('p3')
disp(p3)

p1theory = ones(1, m + 1).' * steadystate(1, 1)
p2theory = ones(1, m + 1).' * steadystate(1, 2)
p3theory = ones(1, m + 1).' * steadystate(1, 3)
p4theory = ones(1, m + 1).' * steadystate(1, 4)
p5theory = ones(1, m + 1).' * steadystate(1, 5)
p6theory = ones(1, m + 1).' * steadystate(1, 6)
p7theory = ones(1, m + 1).' * steadystate(1, 7)
p8theory = ones(1, m + 1).' * steadystate(1, 8)
p9theory = ones(1, m + 1).' * steadystate(1, 9)
p10theory = ones(1, m + 1).' * steadystate(1, 10)
p11theory = ones(1, m + 1).' * steadystate(1, 11)
p12theory = ones(1, m + 1).' * steadystate(1, 12)
p13theory = ones(1, m + 1).' * steadystate(1, 13)
p14theory = ones(1, m + 1).' * steadystate(1, 14)

scf(0); clf;
plot(p1, 'r')
plot(p2, 'g')
plot(p3, 'b')
plot(p4, 'black')
plot(p5, 'yellow')
plot(p1theory, 'r')
plot(p2theory, 'g')
plot(p3theory, 'b')
plot(p4theory, 'black')
plot(p5theory, 'yellow')
title('x[1]:赤, x[2]:緑, x[3]:青, x[4]:黒, x[5]:黄 とsteadystate（横線）の比較');

scf(1);
plot(p6, 'r')
plot(p7, 'g')
plot(p6theory, 'r')
plot(p7theory, 'g')
title('x[6]:赤, x[7]:緑 とsteadystate（横線）の比較');

scf(2);
plot(p8, 'r')
plot(p9, 'g')
plot(p10, 'b')
plot(p8theory, 'r')
plot(p9theory, 'g')
plot(p10theory, 'b')
title('x[8]:赤, x[9]:緑, x[10]:青 とsteadystate（横線）の比較');

scf(3);
plot(p11, 'r')
plot(p12, 'g')
plot(p13, 'b')
plot(p14, 'black')
plot(p11theory, 'r')
plot(p12theory, 'g')
plot(p13theory, 'b')
plot(p14theory, 'black')
title('x[11]:赤, x[12]:緑, x[13]:青, x[14]:黒 とsteadystate（横線）の比較');
