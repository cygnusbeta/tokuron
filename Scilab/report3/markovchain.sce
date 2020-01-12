// Markov 連鎖の定常状態を遷移確率行列 P の固有値・固有ベクトルから求める。
// また、連鎖の差分方程式（漸化式）を反復して、定常状態に収束するか調べる。
clear
xdel(winsid()); // すべてのグラフィックウィンドウを閉じる

// 状態の数と行列 P を与える。
p = [ 0.3 0 0.7; ...
      0 0.6 0.4; ...
      0.1 0.7 0.2]; // テキスト p.43 の例
n = size(p,1);  // 正方行列 p を与えたと仮定している。
disp('n')
disp(n)

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
x = [ 1 0 0 ]; // 初期条件（別のも試してみてください。）
m = 20; // 反復回数
scf(0); clf;

p1 = zeros(1, m + 1)
p2 = zeros(1, m + 1)
p3 = zeros(1, m + 1)

p1(1, 1) = x(1, 1)
p2(1, 1) = x(1, 2)
p3(1, 1) = x(1, 3)

for i=2:m+1
  x = x * p
  steadystate = steadystate / sum(steadystate)
  p1(1, i) = x(1, 1)
  p2(1, i) = x(1, 2)
  p3(1, i) = x(1, 3)
end

disp('p1')
disp(p1)
disp('p2')
disp(p2)
disp('p3')
disp(p3)

plot(p1, 'r')
plot(p2, 'g')
plot(p3, 'b')

title('x[1]:赤, x[2]:緑, x[3]:青 とsteadystate（横線）の比較');

p1theory = ones(1, m + 1).' * steadystate(1, 1)
p2theory = ones(1, m + 1).' * steadystate(1, 2)
p3theory = ones(1, m + 1).' * steadystate(1, 3)
plot(p1theory, 'r')
plot(p2theory, 'g')
plot(p3theory, 'b')
