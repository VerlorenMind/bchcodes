#include "../include/GF.h"

uint64_t mod(uint64_t a, uint64_t b)
{
    /* A proper mod function. */
    while(a >= b)
    {
        a -= b;
    }
    return a;
}

GF2::GF2(uint64_t pow)
{
    /* Generating a log table using primitive polynomial. */
    if(pow > 16 || pow < 2)
    {
        std::cerr<<"Pow for GF is invalid: "<<pow<<std::endl;
    }
    this->deg = pow;
    primitive = PRIMITIVE_POLYS[pow];
    num = uint64_t((1<<pow)-1);

    logs = new std::pair<uint64_t, uint64_t>[num];
    std::pair<uint64_t, uint64_t> gen;
    gen.first = 1;
    gen.second = 0;
    for(uint64_t i=0; i<num; ++i)
    {
        logs[i] = gen;
        ++gen.second;
        gen.first <<= 1;
        if((gen.first >> pow) == 1)
        {
            gen.first ^= this->primitive;
        }
    }
#ifdef DEBUG_GEN
    std::cout<<"Logs table is:\n";
    for(uint64_t i=0; i<num; ++i)
    {
        std::cout<<std::bitset<4>(logs[i].first)<<"\t"<<logs[i].second<<"\n";
    }
#endif
}

GF2::~GF2()
{
    deg = 0;
    primitive = 0;
    logs = nullptr;
    num = 0;
}

uint64_t GF2::add(const uint64_t& a, const uint64_t& b) const
{
    if(a >= (1 << this->deg))
    {
        std::cerr<<"Can't add element: a="<<a<<" is not from GF(2^"<<deg<<")\n";
    }
    if(b >= (1 << this->deg))
    {
        std::cerr<<"Can't add element: b="<<b<<" is not from GF(2^"<<deg<<")\n";
    }
    return a ^ b;
}

uint64_t GF2::mul(const uint64_t& a, const uint64_t& b) const
{
    if(a >= (1 << this->deg))
    {
        std::cerr<<"Can't mul element: a="<<a<<" is not from GF(2^"<<deg<<")\n";
    }
    if(b >= (1 << this->deg))
    {
        std::cerr<<"Can't mul element: b="<<b<<" is not from GF(2^"<<deg<<")\n";
    }
    if (a == 0 || b == 0)
    {
        return 0;
    }
    uint64_t deg1, deg2;
    for(uint64_t i=0; i<num; ++i)
    {
        if(logs[i].first == a)
        {
            deg1 = logs[i].second;
        }
        if(logs[i].first == b)
        {
            deg2 = logs[i].second;
        }
    }
    uint64_t deg = mod(deg1+deg2, num);
#ifdef DEBUG
    std::cout<<"Deg of a: "<<deg1<<"; Deg of b: "<<deg2<<"; Deg of mult:"<<deg<<"; Num+1="<<num+1<<"\n";
#endif
    return logs[deg].first;
}

uint64_t GF2::inv(const uint64_t& a) const
{
    if(a >= (1 << this->deg))
    {
        std::cerr<<"Can't inv element: a="<<a<<" is not from GF(2^"<<deg<<")\n";
    }
    if(a==0)
    {
        return 0;
    }
    uint64_t a_pow;
    for(uint64_t i=0; i<num; ++i)
    {
        if(logs[i].first == a)
        {
            a_pow = logs[i].second;
            break;
        }
    }
    a_pow = num - a_pow;
    return logs[a_pow].first;
}


uint64_t GF2::min_poly(const uint64_t& a) const
{
    if(a >= (1 << this->deg))
    {
        std::cerr<<"Can't find min poly for element: a="<<a<<" is not from GF(2^"<<deg<<")\n";
    }
    /* Finding a minimal polynomial for the element a of GF(2^m) */
    std::vector<uint64_t> pows;
    pows.push_back(a);
    uint64_t atemp = mul(a, a);
    // Calculating powers of a^2 until we get a:
    while(atemp != a)
    {
        pows.push_back(atemp);
        atemp = mul(atemp, atemp);
    }
    std::vector<uint64_t> poly(pows.size()+1);
    std::vector<uint64_t> temp(pows.size()+1);
    poly[0] = pows[0];
    poly[1] = 1;
    // Multiplying (x - bi), where bi - power of a^2
    for(uint64_t i=1; i<pows.size(); ++i)
    {
        for(uint64_t j=0; j<poly.size()-1; ++j)
        {
            temp[j+1] = poly[j];
        }
        for(uint64_t j=0; j<poly.size(); ++j)
        {
            temp[j] = add(temp[j], mul(pows[i], poly[j]));
        }
        poly = temp;
        temp.assign(temp.size(), 0);
    }
    uint64_t ans = 0;
    for(uint64_t i=0; i< poly.size(); ++i)
    {
        ans += (poly[i] << i);
    }
    return ans;
}

GF2::GF2() {
    deg = 0;
    primitive = 0;
    num = 0;
    logs = NULL;
}

bool GF2::is_in(const uint64_t& x) const {
    return x < (1 << deg);
}

uint64_t GF2::get_deg() const {
    return deg;
}

GF2& GF2::operator=(const GF2 &a) {
    this->deg = a.deg;
    this->num = a.num;
    this->primitive = a.primitive;
    delete [] this->logs;
    this->logs = new std::pair<uint64_t, uint64_t>[this->num];
    for(uint64_t i=0; i<this->num; ++i)
    {
        this->logs[i] = a.logs[i];
    }
    return *this;
}

GF2::GF2(const GF2 &a) {
    deg = a.deg;
    num = a.num;
    primitive = a.primitive;
    logs = new std::pair<uint64_t, uint64_t>[this->num];
    for(uint64_t i=0; i<num; ++i)
    {
        logs[i] = a.logs[i];
    }
}

uint64_t GF2::get_elem_deg(uint64_t elem) {
    if(!is_in(elem))
    {
        std::cerr<<"Can't get elem degree; x="<<elem<<" is not from GF(2^"<<deg<<")\n";
    }
    for(uint64_t i=0; i<num; ++i)
    {
        if(logs[i].first == elem)
        {
            return logs[i].second;
        }
    }
    return 0;
}

uint64_t GF2::get_num() const {
    return num;
}


uint64_t mul_poly(uint64_t a, uint64_t b)
{
    uint64_t ans = 0;
    for(uint64_t i=0; b != 0; ++i)
    {
        if(b % 2 == 1)
        {
            ans ^= (a << i);
        }
        b >>= 1;
    }
    return ans;
}

void div_poly(const uint64_t& a, const uint64_t& b, uint64_t &q, uint64_t &r) {
    if (b == 0)
    {
        std::cerr<<"Invalid division, b=0!\n";
        return;
    }
    uint64_t b_pow = 0;
    uint64_t b_temp = b;
    r = a;
    q = 0;
    while(b_temp > 0)
    {
        ++b_pow;
        b_temp >>= 1;
    }
    --b_pow;
    uint64_t r_temp;
    uint64_t r_pow;
    while(b <= r)
    {
        r_pow = 0;
        r_temp = r;
        while(r_temp > 0)
        {
            ++r_pow;
            r_temp >>= 1;
        }
        --r_pow;
        q ^= uint64_t(1) << (r_pow - b_pow);
        r ^= b << (r_pow - b_pow);
    }
}

uint64_t gcd_poly(uint64_t a, uint64_t b) {
    uint64_t q, r;
    while(b != 0)
    {
        div_poly(a, b, q, r);
        a = b;
        b = r;
    }
    return a;
}

uint64_t lcm_poly(uint64_t a, uint64_t b)
{
    uint64_t denom = gcd_poly(a, b);
    uint64_t ab = mul_poly(a, b);
    uint64_t q, r;
    div_poly(ab, denom, q, r);
    return q;
}

uint64_t deg_poly(uint64_t poly)
{
    uint64_t ans = 0;
    while(poly > 0)
    {
        ++ans;
        poly >>= 1;
    }
    --ans;
    return ans;
}

GF2X::GF2X(): deg(1), F(2) {
    c = new uint64_t[1];
    c[0] = 0;
}

GF2X::GF2X(const GF2& field, uint64_t *coeffs, const uint64_t& pow) {
    F = field;
    this->deg = pow;
    for(uint64_t i=0; i<=deg; ++i)
    {
        if(!F.is_in(coeffs[i]))
        {
            std::cerr<<"Invalid coeffs for GF2X: "<<i<<"th coeff="<<coeffs[i]<<" is not from GF(2^"<<F.get_deg()<<")\n";
        }
        break;
    }
    c = coeffs;
}

GF2X::GF2X(const GF2& field, const std::vector<uint64_t>& coeffs) {
    F = field;
    deg = uint64_t(coeffs.size()-1);
    c = new uint64_t[deg+1];
    for(uint64_t i=0; i<=deg; ++i)
    {
        if(!F.is_in(coeffs[i]))
        {
            std::cerr<<"Invalid coeffs for GF2X: "<<i<<"th coeff="<<coeffs[i]<<" is not from GF(2^"<<F.get_deg()<<")\n";
        }
        break;
    }
    for(uint64_t i=0; i<=deg; ++i)
    {
        c[i] = coeffs[i];
    }
}

GF2X::~GF2X() {
    delete [] c;
}

GF2X GF2X::operator+(const GF2X &a) const {
    if(this->F.get_deg() != a.F.get_deg())
    {
        std::cerr<<"Can't sum two polys: fields does not match\n";
        return GF2X();
    }
    uint64_t* ans;
    uint64_t max_pow = a.deg > this->deg ? a.deg : this->deg;
    ans = new uint64_t[max_pow+1];
    for(uint64_t i=0; i<=max_pow; ++i)
    {
        if(i <= a.deg && i <= this->deg)
        {
            ans[i] = F.add(a.c[i], this->c[i]);
        }
        else if(i > a.deg)
        {
            ans[i] = this->c[i];
        }
        else if(i > this->deg)
        {
            ans[i] = a.c[i];
        }
    }
    GF2X ans_p(F, ans, max_pow);
    ans_p.cut();
    return ans_p;
}

GF2X GF2X::operator*(const GF2X &a) const {
    if(this->F.get_deg() != a.F.get_deg())
    {
        std::cerr<<"Can't mul two polys: fields does not match\n";
        return GF2X();
    }
    uint64_t ans_pow = this->deg + a.deg;
    uint64_t* ans = new uint64_t[ans_pow+1];
    for(uint64_t i=0; i<=ans_pow; ++i)
    {
        ans[i] = 0;
    }
    for(uint64_t i=0; i<=a.deg; ++i)
    {
        for(uint64_t j=0; j<=this->deg; ++j)
        {
            ans[i+j] = F.add(ans[i+j], F.mul(a.c[i], this->c[j]));
        }
    }
    GF2X ans_p(F, ans, ans_pow);
    ans_p.cut();
    return ans_p;
}

uint64_t GF2X::eval(const uint64_t &x) const {
    if(!F.is_in(x))
    {
        std::cerr<<"Can't eval poly: x="<<x<<" is not from GF(2^"<<F.get_deg()<<")\n";
        return 0;
    }
    uint64_t ans = 0;
    for(int i = this->deg; i >= 0; --i)
    {
        ans = F.mul(x, ans);
        ans = F.add(ans, c[i]);
    }
    return ans;
}

std::string GF2X::print() {
    std::string ans = "";
    std::string x = "";
    for(int i=0; i<=deg; ++i)
    {
        if(c[i] == 0)
        {
            continue;
        }
        else if(c[i] == 1)
        {
            x = "x^" + std::to_string(i);
        }
        else
        {
            x = std::to_string(c[i]) + "x^" + std::to_string(i);
        }
        ans = ans == "" ? x : x + " + " + ans;
    }
    return ans;
}

uint64_t &GF2X::operator[](uint64_t i) {
    return *(c+i);
}

void GF2X::cut() {
    uint64_t old_deg = deg;
    while(c[deg] == 0 && deg > 0)
    {
        --deg;
    }
    if(old_deg != deg)
    {
        uint64_t* ctemp = new uint64_t[deg+1];
        for(uint64_t i=0; i<=deg; ++i)
        {
            ctemp[i] = c[i];
        }
        delete [] c;
        c = ctemp;
    }
}

GF2X::GF2X(const GF2X &x) {
    deg = x.deg;
    F = x.F;
    c = new uint64_t[deg+1];
    for(uint64_t i=0; i<=deg; ++i)
    {
        c[i] = x.c[i];
    }
}

GF2X &GF2X::operator=(const GF2X &x) {
    deg = x.deg;
    F = x.F;
    if(c != NULL) delete [] c;
    c = new uint64_t[deg+1];
    for(uint64_t i=0; i<=deg; ++i)
    {
        c[i] = x.c[i];
    }
    return *this;
}

GF2X::GF2X(const GF2 &field, uint64_t poly) {
    F = field;
    deg = 0;
    if (poly > 0) {
        uint64_t temp = poly;
        while (temp > 0) {
            ++deg;
            temp >>= 1;
        }
        --deg;
    }
    c = new uint64_t[deg+1];
    for(uint64_t i=0; i<=deg; ++i)
    {
        c[i] = poly % 2;
        poly >>= 1;
    }
}

GF2X GF2X::operator*(uint64_t x) const {
    GF2X answer(*this);
    for(uint64_t i=0; i<=deg; ++i)
    {
        answer.c[i] = F.mul(c[i], x);
    }
    answer.cut();
    return answer;
}

uint64_t GF2X::get_deg() {
    return deg;
}

uint64_t* GF2X::roots(uint64_t& num)
{
    uint64_t alpha = 1;
    uint64_t* root = new uint64_t(F.get_deg());
    num=0;
    for(uint64_t j=0; j<F.get_num(); ++j)
    {
        if(eval(alpha) == 0)
        {
            root[num] = alpha;
            ++num;
        }
        alpha = F.mul(alpha, 2);
    }
    if(num < F.get_deg())
    {
        uint64_t* temp = new uint64_t[num];
        for(uint64_t i=0; i<num; ++i)
        {
            temp[i] = root[i];
        }
        delete [] root;
        root = temp;
    }
    return root;
}


