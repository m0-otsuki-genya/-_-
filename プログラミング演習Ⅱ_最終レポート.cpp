#include<string>
#include<stdio.h>
#include<vector> 
#include<iostream>
#include<cmath>
#include<algorithm>
#include<list>
#include<numeric>
using namespace std;

double I1(double& x, double& y) {
	cout <<"断面二次モーメントは、" << (x * x * x * y) / 12 << "[m^4]\n";
	return (x * x * x * y) / 12;
}
double I2(double& x, double& y) {
	cout << "断面二次モーメントは、" << ((x * x * x * x) - (y * y * y * y)) * 3.14 / 64 << "[m^4]\n";
	return ((x * x * x * x) - (y * y * y * y)) * 3.14 / 64;
}
double BMmax(double& x, double& y, double& z) {
	return (x* y) / (2 * z);
}

int main() {
	vector<double> vec(6);
	
	cout <<"片持ちばり？それとも両端支持ばり？前者は1を、後者は2を押して。\n";
	string str;
	cin >> str;
	if (str == "1") {
		cout << "片持ちばりね？";
		vec[0] = 1;
	}
	else if (str == "2") {
		cout << "両端支持ばりね？";
		vec[0] = 2;
	}
	else cout << "1と2以外を入力したね？君。最初からやり直してくれ。～完～";

	cout << "断面形状について尋ねるよ。長方形（正方形）断面なら1を、円形断面（中実、中空）なら2を押して。\n";
	cin >> str;
	if (str == "1") {
		cout << "長方形（正方形）断面ね。";
		vec[1] = 1;
	}
	else if (str == "2") {
		cout << "円形断面ね？";
		vec[1] = 2;
	}
	else cout << "1と2以外を入力したね？君。最初からやり直してくれ。～完～";

	if (vec[1] == 1) {
		cout << "長方形断面の幅bと高さhを聞くよ。まずh高さをm単位で教えてね。\n";
		string str;
		cin >> str;
		double a = stod(str);
		vec[2] = a;
		cout << "続いて幅bをm単位で教えてね。\n";
		cin >> str;
		double b = stod(str);
		vec[3] = b;
	}
	if (vec[1] == 2) {
		cout << "円形断面の外径と内径を聞くよ。まず外径をm単位で教えてね。\n";
		cin >> str;
		double a = stod(str);
		vec[2] = a;
		cout << "続いて内径をm単位で教えてね。中実丸棒の場合0と答えてな。\n";
		cin >> str;
		double b = stod(str);
		vec[3] = b;
//		if (a <= b) cout<<"外径より内径のが大きいぞ、どういうこっちゃ！\n" << endl;
	}

	vector<double> vecx;
	std::vector<double> vecf;
	vector<double> vecdx;
	vector<double> vecdf;
	vector<double> vectx;
	vector<double> vectf;
	vector<double> vec0xl;
	vec0xl.push_back(0);

	do {
		cout << "荷重について聞くよ。集中荷重を追加するなら1を、等分布荷重を追加するなら2を、三角形の分布荷重は3を押して。\n";
		cin >> str;
		if (str == "1") {
			cout << "集中荷重ね。では荷重[N]とその位置x[m]を聞こう。まず荷重[N]を教えて。　";
			cin >> str;
			double c = stod(str);
			vecf.push_back(c);
			cout << "\n続いて位置x[m]を教えて。(固定端をlとして、x=0～l)　";
			cin >> str;
			double d = stod(str);
			vecx.push_back(d);
			vec0xl.push_back(d);
		}
		else if (str == "2") {
			cout << "等分布荷重ね。どこからどこまで、何荷重[N/m]か聞く。まずは開始点[m]を教えて。";
			cin >> str;
			double c = stod(str);
			vecdx.push_back(c);
			vec0xl.push_back(c);
			cout << "次に荷重の終点の座標[m]を教えて。";
			cin >> str;
			double d = stod(str);
			vecdx.push_back(d);
			vec0xl.push_back(d);
			cout << "最後に等荷重[N/m]を教えて。";
			cin >> str;
			double f = stod(str);
			vecdf.push_back(f);
		}
		else if (str == "3") {
			cout << "三角形の分布荷重ね。じゃあどこからどこまで、何荷重[N/m]か聞く。まずは一番xの小さい点[m]を教えて。";
			cin >> str;
			double c = stod(str);
			vec0xl.push_back(c);
			cout << "次に,三角形の頂点の座標[m]を教えて。";
			cin >> str;
			double d = stod(str);
			vec0xl.push_back(d);
			cout << "最後に、頂点での三角形の高さ[N/m]を教えて。";
			cin >> str;
			double g = stod(str);
			vectf.push_back(g);
		}
		else cout << "1と2以外を入力したね？君。最初からやり直してくれ。～完～";

		cout << "\n他の荷重を追加する？するならaddと、しないならそれ以外の何かを適当に入力して。\n";
		cin >> str;
	} while (str == "add");

	cout << "最後に、はりの全長[m]を教えて。\n";
	cin >> str;
	double h = stod(str);
	vec[4] = h;
	vec0xl.push_back(h);
	cout << "\n";


	vector<double> vecDf;
	vector<double> vecTf;
	for (int i = 0; i < vecdf.size();i++) {
		vecDf.push_back(vecdf[i]*(vecdx[1+2*i] - vecdx[2*i]));
	}
	for (int i = 0; i < vectf.size(); i++) {
		vecTf.push_back(vectf[i] * (vecdx[1 + 2 * i] - vecdx[2 * i])/2);
	}


	double Rb = 0;
	for (int i = 0; i < vecf.size(); i++) {
		for (int j = 0; j < vecDf.size(); j++) {
			for (int k = 0; k < vecTf.size(); k++) {
				Rb += ((vecf[i] * vecx[i]) + (vecDf[j] * (vecdx[2 * j] + vecdx[1 + 2 * j]) / 2) + (vecTf[k] * (vectx[2 * k] + 2 * vectx[1 + 2 * k]) / 3)) / vec[4];
			}
		}
	}
	if (vec[0] == 2) {
		double Ra = accumulate(vecf.begin(), vecf.end(),0) + accumulate(vecDf.begin(), vecDf.end(),0) + accumulate(vecTf.begin(), vecTf.end(),0) - Rb;
	}
	else {
		double Ra = 0;
	}


	double Mx = 0;
	int xmax=0;
	auto F1max = max_element(vecf.begin(), vecf.end());
	auto F2max = max_element(vecDf.begin(), vecDf.end());
	auto F3max = max_element(vecTf.begin(), vecTf.end());
	if (F1max > F2max && F1max > F3max) {
		int ia = std::distance(vecf.begin(),F1max);
		xmax = ia;
		vec.push_back(vecx[ia]);
			for (int i=0; i<xmax; i++) {
			Mx += vecf[i]*vecx[i];
		}
	}
	else if (F2max > F1max && F2max > F3max) {
		int ib = std::distance(vecDf.begin(), F2max);
		xmax = ib;
		vec.push_back((vecdx[2*ib]+vecdx[1+2*ib])/2);
		if (xmax % 2 == 0) {
			for (int i=0; i < xmax; i++) {
				Mx += -vecDf[i] * (vecdx[2*i]+vecdx[1+2*i])/2;
			}
	}
		else if (xmax % 2 == 1) {
			Mx += -vecDf[xmax] * (vecdx[2*xmax] + vec[5]) / 2;
			for (int i=0; i < (xmax - 1); i++) {
				Mx += -vecDf[i] * (vecdx[2 * i] + vecdx[1 + 2 * i]) / 2;
			}
		}
	}
	else if (F3max > F1max && F3max > F2max) {
		int ic = std::distance(vecTf.begin(), F3max);
		xmax = ic;
		vec.push_back((vectx[2*ic]+2*vectx[1+2*ic])/3);
		if (xmax % 2 == 0) {
			for (int i=0; i < xmax; i++) {
				Mx += -vecTf[i] * (vectx[2 * i] + 2*vectx[1 + 2 * i]) / 3;
			}
		}
		else if (xmax % 2 == 1) {
			Mx += -vecTf[xmax] * (vectx[2 * xmax] + 2*vec[5]) / 3;
				for (int i=0; i < (xmax - 1); i++) {
					Mx += -vecTf[i] * (vectx[2 * i] + 2*vectx[1 + 2 * i]) / 3;
				}
		}
	}

	if (vec[1] == 1) vec.push_back(I1(vec[2], vec[3]));
	else if(vec[1] == 2) vec.push_back(I2(vec[2], vec[3]));


	if (vec[0] == 1) {
		double Ml = 0;
		for (int i = 0; i < vecf.size(); i++) {
			for (int j = 0; j < vecDf.size(); j++) {
				for (int k = 0; k < vecTf.size(); k++) {
					Ml += -vecf[i]*(vec[4]-vecx[i]) - vecDf[j]*(vec[4]-vecdx[2*j]-(vecdx[1+2*j]-vecdx[2*j]) / 2) - vecTf[k]*(vec[4]-vectx[2*k]-2*(vectx[1+2*k]-vectx[2*k]) / 3);
				}
			}
		}
		BMmax(Ml,vec[2], vec[6]);
		cout << "x=" << vec[4] << "[m](固定端)で曲げ応力は最大となる。\n";
		cout << "はりに生ずる最大曲げ応力は、" << BMmax << "[Pa]である。";
	}
	else if (vec[0] == 2) {
		BMmax(Mx, vec[2], vec[6]);
		cout << "x=" << vec[5] << "[m]で曲げ応力は最大となる。\n";
		cout << "はりに生ずる最大曲げ応力は、" << BMmax << "[Pa]である。";
	}


	for (auto i = vec.begin(); i != vec.end(); ++i) {
		cout << *i << ",　";
	}
	cout << endl;
}

// プログラムの実行: Ctrl + F5 または [デバッグ] > [デバッグなしで開始] メニュー
// プログラムのデバッグ: F5 または [デバッグ] > [デバッグの開始] メニュー