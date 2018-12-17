#include<iostream>
#include<functional>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<vector>
#include<tuple>
#include<random>
#include<algorithm>
#include<regex>
#include<bitset>
#include<thread>
#include<string>
#include<chrono>
#include<unordered_map>
#include<stack>
#include<boost/range/irange.hpp>
#include<boost/range/join.hpp>
#include"magicmoves.cpp"
#include"magicmoves.hpp"
#include <strings.h>
#include<array>
#include<mutex>
#include<future>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>
int ffs(int i);
#define _GNU_SOURCE
#include <string.h>
int ffsl(long int i);
int ffsll(long long int i);

using namespace std;

template<typename T> int f_set(T t){
  return ffsll(t) - 1;
}

namespace util{

	template<typename R> R rand_real(){
	  return R(rand() % (1 << 16)) / R(1 << 16);
	}

  template<typename I> struct facade{
    using iterator = I;
    I first;
    I last;
    I begin(){
      return first;
    }
    I end(){
      return last;
    }
    facade(I&& f, I&& l) : first(forward<I>(f)), last(forward<I>(l)){}
  };
  
  template<typename It, typename Trans> auto best_element(It&& a, It&& b, Trans&& func){
    auto best = make_tuple(a, func(*a));
    size_t count(0);
    for(; a != b; ++a){
      if(auto temp = func(*a); temp >= get<1>(best)){
        if(temp == get<1>(best)){
          ++count;
          if((rand() % count) == 0) best = make_tuple(a, temp);
        }else best = make_tuple(a, temp);
      }
    }
    return best;
  }
  
  template<typename I> facade<I> make_facade(I a, I b){
    return facade<I>(forward<I>(a), forward<I>(b));
  }
  
  struct bit_iterator{
    using U64 = unsigned long long int;
    U64 local;
    size_t operator*(){
      return f_set(local);
    }
    
    bool operator==(const bit_iterator& other){
      return other.local == local;
    }
    
    bool operator!=(const bit_iterator& other){
      return !(*this == other);
    }
    
    bit_iterator& operator++(){
      local ^= (U64(1) << f_set(local));
      return *this;
    }    
    bit_iterator(U64 in) : local(in) {}
  };
  
  template<typename T> auto serialize(T&& t){
    return make_facade(bit_iterator(forward<T>(t)), bit_iterator(0));
  }
  
  template<size_t N> void rotate(bitset<N>& input){
    for(size_t i(0); i < N / 2; ++i){
      bool temp = input[N - i - 1];
      input[N - i - 1] = input[i];
      input[i] = temp;
    }
  }
  
  template<typename T> void rotate(T& input){
    bitset<64> temp(input);
    rotate(temp);
    input = temp.to_ulong();
  }
  
  double sigmoid(double x){
    return 2.0 / (1.0 + pow(2, -x));
  }
 
  template<typename T> auto ensure(size_t count, T f){
    return [=](auto a, auto b){
      size_t neg(0);
      size_t pos(0);
      for(size_t i(0); i < count; ++i){
        if(f(a, b)) ++pos;
        else ++neg;
      }
      for(size_t i(0); i < count; ++i){
        if(f(b, a)) ++neg;
        else ++pos;
      }
      return (pos > neg);
    };
  }
  
}

namespace cartes{

	template<typename R> R rand_real(){
	  return R(rand() % (1 << 16)) / R(1 << 16);
	}

  template<typename I, typename T> struct node{
    I _id = 0;
    I _1 = 0;
    I _2 = 0;
    optional<T> result{};
    
    string to_string(){
      string temp{};
      temp = std::to_string(_id) + " _ " + std::to_string(_1) + " " + std::to_string(_2);
      return temp;        
    };
    
  };

  template<size_t L, size_t O, size_t C, typename T> struct program{
    using value_type = T;
    static const size_t max_index = L + L + O;
    
    vector<function<T(T, T)>> grammar;
    
    
    array<T, O> input;
    array<T, C> constants;
    //vector enables for larger network at the cost of reduced performance
    vector<node<size_t, T>> genome = vector<node<size_t, T>>(L);
    
    T index(size_t l){
      if(l < O) return input[l];
      if(l < (O + C)) return constants[l - O];
      else return exec(genome[l - C - O]);
    }
    
    program<L, O, C, T>& wipe(){
      for(auto& elem : genome) elem.result = {};
      return *this;
    }
    
    T operator()(array<T, O> i){
      input = i;
      T temp = exec(genome.back()); wipe();
      return temp;
    }
    
    T exec(node<size_t, T>& input){
      if(input.result.has_value()) return input.result.value(); 
      else{
        input.result = grammar[input._id](index(input._1), index(input._2));
        return input.result.value();
      }
    }
    
    string to_string(){
      string result;
      for(auto& n : genome) result += n.to_string() + '\n';
      result += "\n\n";
      for(auto& c : constants) result += c.to_string() + '\n';
      return result;
    }
    
    template<typename R> program<L, O, C, T>& permute(R pr){
      size_t l(0);
      for(node<size_t, T>& n : genome){
        if(rand_real<R>() < pr) n._id = rand() % grammar.size();
        //?? 1/6 probability of solely referencing input parameter ... questionable
        if(rand_real<R>() < pr) n._1 = ((rand() % 18) == 0) ? (rand() % (l + C + O)) : (rand() % O);
        if(rand_real<R>() < pr) n._2 = ((rand() % 18) == 0) ? (rand() % (l + C + O)) : (rand() % O);
        ++l;
      }
      for(T& elem : constants) elem.permute(forward<R>(pr));
      return *this;
    }
    
    static program<L, O, C, T> create(vector<function<T(T, T)>> grammar){
      return program<L, O, C, T>(grammar).permute(1.0);
    }
    
    template<typename I> static program<L, O, C, T> load(vector<function<T(T, T)>> grammar, I& input){
      program<L, O, C, T> result(grammar);
      for(auto& node : result.genome){
        string delim;
        input >> node._id;
        input >> delim;
        input >> node._1;
        input >> node._2;
      }
      for(T& constant : result.constants) constant.load(input);
      return result;
    }
    
    program(vector<function<T(T, T)>> gram) : grammar(gram){}
    program(){}
  };

}

namespace ga{

  template<typename T, size_t P> struct population{
    array<T, P> pop;
    T seed;
    
    template<typename Cmp> population<T, P> find_seed(Cmp&& GT){
      seed = pop.back();
      for(T& prog : pop){
        if(GT(prog, seed)) seed = prog;
      }
      return *this;
    }
    
    template<typename Cmp> population<T, P>& tournament_seed(Cmp&& GT, size_t times = 1){
      array<int, P> scores; scores.fill(0);
      
      for(size_t i(0); i < P; ++i){
        for(size_t j(0); j < P; ++j){
          for(size_t k(0); k < times; ++k){
            if(i != j){
              auto temp = GT(pop[i], pop[j]);
              scores[i] += temp;
              scores[j] -= temp;
            }
          }
        }
      }
      for(auto& score : scores) cout << score << " ,";
      seed = pop[distance(scores.begin(), max_element(scores.begin(), scores.end()))];
      return *this;
    }
    
    template<typename Cmp> population<T, P>& challenge_seed(Cmp&& GT, size_t times = 1){
      auto score = [&, this](T& in){
        int val(0);
        bool side = true;
        for(size_t i(0); i < times; ++i, side = !side){
          if(side) val += GT(in, seed);
          else val -= GT(seed, in);
        }
        return val;
      };
      optional<T> best; int result(0);
      for(auto& prog : pop){
        if(!best.has_value()){
          best = prog; result = score(prog);
        }else if(int temp = score(prog); temp > result){
            best = prog; result = temp;
        }
      }
      cout << "best score :: " << result << endl;
      seed = best.value();
      return *this;
    }
    
    template<typename R> population<T, P>& evolve(R pr, R replacement = R(0)){
      for(T& prog : pop){
        if(util::rand_real<R>() > replacement){
          prog = seed;
          prog.permute(pr);
        }
      }
      return *this;
    }
    
    template<typename G> population(G g){
      for(T& prog : pop) prog = T::create(g);
      seed = pop.back();
    }
  };

}

namespace chess{
  //this move generator is solely for training as it is has poor performance,
  //though it probably won't prove too much of a bottleneck anyways. Amadahl's law...

  enum piece_t{KING, QUEEN, ROOK, BISHOP, KNIGHT, PAWN, ALL, EMPTY};
  const array<piece_t, 6> pieces{KING, QUEEN, ROOK, BISHOP, KNIGHT, PAWN};  
  using move = struct{
    piece_t type;
    piece_t capture;
    size_t from;
    size_t to;
  };
  using U64 = unsigned long long int;
  
  struct sq_t{
    int r;
    int f;
    
    template<typename T> sq_t& add(T&& t){
      auto [a, b] = t;
      r += a; f += b;
      return *this;
    }
    
    sq_t(int a, int b) : r(a), f(b){}
  };
  
  struct ray_t{
    int r;
    int f;
    ray_t(int a, int b) : r(a), f(b){}
  };
  
  sq_t inflate(int input){
    return sq_t(input / 8, input % 8);
  }
  
  int deflate(sq_t s){
    return s.r * 8 + s.f;
  }

  bool s_is_valid(sq_t s){
    return (s.r >= 0 && s.f >= 0 && s.r < 8 && s.f < 8);
  }
  
  bool is_valid(int s){
    return 0 <= s && s < 64;
  }
  
  string name(sq_t s){
    array<string, 8> files{{"a", "b", "c", "d", "e", "f", "g", "h"}};
    reverse(files.begin(), files.end());
    reverse(files.begin(), files.end()); //?
    return files[s.f] + to_string(s.r + 1);
  }
  
  string name(int s){
    return name(inflate(s));
  }
  
  string name(move m){
    return name(m.from) + name(m.to);
  }
  

  inline constexpr U64 b_file(int i){
    return U64(72340172838076673) << i;
  }
  
  inline constexpr U64 b_rank(int i){
    return U64(255) << (i * 8);
  }

  
  constexpr U64 black_OOO_rook(){
    return U64(1) << U64(63);
  }
  
  constexpr U64 white_OOO_rook(){
    return U64(1) << U64(7);
  }
  
  constexpr U64 white_OO_rook(){
    return U64(1);
  }
  
  constexpr U64 black_OO_rook(){
    return U64(1) << U64(56);
  }
  
  template<size_t N> auto offset_tbl(array<ray_t, N> offsets){  
    array<U64, 64> result;
    for(size_t i(0); i < 64; ++i){
      for(auto elem : offsets){
        if(s_is_valid(inflate(i).add(elem))){
          result[i] |= (U64(1) << deflate(inflate(i).add(elem)));
        }
      }
    }
    return result;
  }
  
  namespace knight{
    const array<U64, 64> tb = offset_tbl(array<ray_t, 8>{{
      ray_t(1, 2),
      ray_t(-1, 2),
      ray_t(2, -1),
      ray_t(-2, -1),
      ray_t(-2, 1),
      ray_t(1, -2),
      ray_t(2, 1),
      ray_t(-1, -2)     
    }});  
  }
  
  namespace king{
    const array<U64, 64> tb = offset_tbl(array<ray_t, 8>{{
      ray_t(1, 0),
      ray_t(-1, 0),
      ray_t(0, -1),
      ray_t(1, -1),
      ray_t(1, 1),
      ray_t(0, 1),
      ray_t(-1, -1),
      ray_t(-1, 1)     
    }});
  }

	namespace black_pawn{
		const array<U64, 64> advance = offset_tbl(array<ray_t, 1>{{
			ray_t(-1, 0)
		}});

		const array<U64, 64> attack = offset_tbl(array<ray_t, 2>{{
			ray_t(-1, 1),
			ray_t(-1, -1),
		}});
	}

	namespace white_pawn{
		const array<U64, 64> advance = offset_tbl(array<ray_t, 1>{{
			ray_t(1, 0)
		}});

		const array<U64, 64> attack = offset_tbl(array<ray_t, 2>{{
			ray_t(1, 1),
			ray_t(1, -1),
		}});
	}
	
  
  void print_as_board(U64 in){
    for(int i = 0; i < 8; ++i){
      for(int j = 0; j < 8; ++j){
        cout << (in & (U64(1) << (i * 8 + j)) ? "1 " : ". ");
      }
      cout << '\n';
    }
    cout << '\n';
  }
  
  struct state{
    U64 k = 0;
    U64 q = 0;
    U64 r = 0;
    U64 b = 0;
    U64 n = 0;
    U64 p = 0;
    U64 all = 0;
    U64 enpassant = 0;
		U64 lastdouble = 0;
    bool OOO = true;
    bool OO = true;
    
    U64& operator[](size_t in){
      switch(in){
        case KING: return k;
        case QUEEN: return q;
        case ROOK: return r;
        case BISHOP: return b;
        case KNIGHT: return n;
        case PAWN: return p;
        case ALL: return all;
      }
    }   
    
    template<typename F> constexpr state& apply(F&& f){
      f(k); f(q); f(r); f(b); f(n); f(p); f(all);
      return *this;
    }
    
    state& update(){
      all = 0;
      return apply([this](U64 i){ all |= i; });
    }
    
    static constexpr state black(){
      auto ret = white();
      ret.apply(util::rotate<U64>);
      auto temp = ret.q;
      ret.q = ret.k; ret.k = temp;
      return ret;
    }
    
    static constexpr state white(){
      state ret;
      ret.all = b_rank(1) | b_rank(0);
      ret.p = b_rank(1);
      ret.r = (U64(1) << 7) | 1;
      ret.n = (U64(1) << 1) | (U64(1) << 6);
      ret.b = (U64(1) << 2) | (U64(1) << 5);
      ret.k = U64(1) << 3;
      ret.q = U64(1) << 4;
      return ret;    
    }
  };
  
  bool operator==(state& a, state& b){
    if(a.all != b.all) return false;
    for(piece_t p : pieces){
      if(a[p] != b[p]) return false;
    }
    if(a.enpassant != b.enpassant) return false;
    if(a.lastdouble != b.lastdouble) return false;
    if(a.OO != b.OO) return false;
    if(a.OOO != b.OOO) return false;
    return true; 
  }
  
  struct special_move{
    U64 eliminated = 0;
    state result;
    special_move(state input) : result(input){}
    special_move(U64 el, state input) : eliminated(el), result(input){}
  };
  
  struct board{
    state white = state::white();
    state black = state::black();
    bool white_pov{true};
    size_t count{0};

    state us_start(){ return white_pov ? state::white() : state::black(); }

		U64 b_moves(size_t index){
			return Bmagic(index, white.all | black.all) & ~white.all;
		}
		
		U64 r_moves(size_t index){
			return Rmagic(index, white.all | black.all) & ~white.all;
		}

		U64 q_moves(size_t index){
			return r_moves(index) | b_moves(index);
		}

		U64 n_moves(size_t index){
			return knight::tb[index] & ~white.all;
		}
		
		U64 k_moves(size_t index){
			return king::tb[index] & ~white.all;
		}

		U64 p_advance(size_t index){
			if(white_pov) return white_pawn::advance[index] & ~(black.all | white.all);
			else return black_pawn::advance[index] & ~(black.all | white.all);
		}

		U64 p_attack(size_t index){
			if(white_pov) return white_pawn::attack[index] & black.all;
			else return black_pawn::attack[index] & black.all;
		}

		U64 p_moves(size_t index){
			return p_attack(index) | p_advance(index);
		}

    U64 standard_moves(piece_t p, size_t index){
      switch(p){
        case KING: return k_moves(index);
        case QUEEN: return q_moves(index);
        case ROOK: return r_moves(index);
        case BISHOP: return b_moves(index);
        case KNIGHT: return n_moves(index);
        case PAWN: return p_moves(index);
      }
    }

    bool attacked(size_t index){
      for(piece_t p : pieces){
        if(standard_moves(p, index) & black[p]) return true;
      }
      return false;
    }
    
    bool is_check(){
      bool result = attacked(log2(white.k));
      return result;
    }

    template<typename F> board& double_moves(F&& visitor){
      U64 pawn_rank = white_pov ? b_rank(1) : b_rank(6);
      U64 empty = white_pov ? 
        ((U64(1) << U64(24)) | U64(1) << U64(16)):
        ((U64(1) << U64(32)) | U64(1) << U64(40));
      int offset = (white_pov ? 16 : -16);
      int mask_offset = white_pov ? 8 : -8;
      U64 occ = white.all | black.all;
      for(int from : util::serialize(white.p & pawn_rank)){
        if(!bool((empty << U64(from % 8)) & occ)){
          state cpy = white;
          cpy.enpassant = U64(1) << U64(mask_offset + from);
					cpy.lastdouble = U64(1) << U64(offset + from);
          cpy.p &= ~(U64(1) << U64(from)); cpy.p |= (U64(1) << U64(from + offset));
          cpy.all &= ~(U64(1) << U64(from)); cpy.all |= (U64(1) << U64(from + offset));          
          visitor(special_move(cpy));
        }
      }
      return *this;
    }
    
	  template<typename F> board& enpassant_moves(F&& visitor){
	  	for(size_t from : util::serialize(white.p)){
	  		U64 enattack = (white_pov ? white_pawn::attack : black_pawn::attack)[from] & black.enpassant;
	  		if(bool(enattack)){
	  			state cpy = white;
	  			cpy.p &= ~(U64(1) << from);
	  			cpy.p |= black.enpassant;
	  			cpy.update();
	  			visitor(special_move(black.lastdouble, cpy));
	  		}
	  	}
	  	return *this;
	  }

    template<typename F> board& castle_moves(F&& visitor){
      if(is_check()) return *this;
      auto clear = [this](auto index_range){
        U64 occ = white.all | black.all;
        for(size_t val : index_range){
          if(attacked(val) | (occ & (U64(1) << U64(val)))) return false;
        }
        return true;
      };
      auto context_OOO = white_pov ? boost::irange(4, 7) : boost::irange(60, 63);
      auto context_OO = white_pov ? boost::irange(1, 3) : boost::irange(57, 59);
      if(white.OO && clear(context_OO)){
        state cpy = white;
        if(white_pov){
          cpy.k = U64(1) << U64(1);
          cpy.r &= ~U64(1);
          cpy.r |= U64(1) << U64(2);
        }else{
          cpy.k = U64(1) << U64(57);
          cpy.r &= ~(U64(1) << U64(56));
          cpy.r |= U64(1) << U64(58);
        }
        cpy.update();
        cpy.OO = cpy.OOO = false;
        visitor(special_move(cpy));
      }
      if(white.OOO && clear(context_OOO)){
        state cpy = white;
        if(white_pov){
          cpy.k = U64(1) << U64(5);
          cpy.r &= ~(U64(1) << U64(7));
          cpy.r |= U64(1) << U64(4);
        }else{
          cpy.k = U64(1) << U64(61);
          cpy.r &= ~(U64(1) << U64(63));
          cpy.r |= U64(1) << U64(60);
        }
        cpy.OO = cpy.OOO = false;
        cpy.update();
        visitor(special_move(cpy));
      }
      return *this;
    }
    
    template<typename F> board& promotion_moves(F&& visitor){
      U64 final_rank = white_pov ? b_rank(6) : b_rank(1);
      for(size_t from : util::serialize(white.p & final_rank)){
        for(size_t to : util::serialize(p_moves(from))){
          if((to / 8 == 0) || (to / 8 == 7)){
              
              auto type = QUEEN;
              if(type != PAWN && type != KING){
                state cpy = white;
                cpy.p &= ~(U64(1) << from);
                cpy.all &= ~(U64(1) << from);
                cpy.all |= (U64(1) << to);
                cpy[type] |= (U64(1) << to);
                visitor(special_move(cpy));
              }

          }
        }
      }
      return *this;
    }

		template<typename F> board& conv_moves(F&& visitor){
      U64 final_rank = white_pov ? b_rank(6) : b_rank(1);
      auto helper = [&visitor](piece_t p, const U64& b, auto generator){
        for(size_t i : util::serialize(b)){
          for(size_t j : util::serialize(generator(i))) visitor(move{p, EMPTY, i, j});
        }
      };
      helper(KING, white.k, [this](auto&& in){return k_moves(in);});
      helper(QUEEN, white.q, [this](auto&& in){return q_moves(in);});
      helper(ROOK, white.r, [this](auto&& in){return r_moves(in);});
      helper(BISHOP, white.b, [this](auto&& in){return b_moves(in);});
      helper(KNIGHT, white.n, [this](auto&& in){return n_moves(in);});
      helper(PAWN, white.p, [this](auto&& in){return p_moves(in) & ~(b_rank(0) | b_rank(7));});
      
      return *this;
		}

    template<typename F> board& visit_pseudo_moves(F&& visitor){
      return conv_moves(forward<F>(visitor)).double_moves(forward<F>(visitor))
             .castle_moves(forward<F>(visitor)).promotion_moves(forward<F>(visitor)).enpassant_moves(forward<F>(visitor));
    }
    
    template<typename F> board& visit_positions(F&& f){
      return visit_pseudo_moves([&f, this](auto&& act){
        if(auto temp = do_(act); !temp.is_check()) f(temp.flip());
      });
    }

    board do_(const move& m){
      auto context_OOO = [this](bool pov){ return pov ? white_OOO_rook() : black_OOO_rook(); };
      auto context_OO = [this](bool pov){ return pov ? white_OO_rook() : black_OO_rook(); };   
       
      U64 from = (U64(1) << m.from);
      U64 to = (U64(1) << m.to);
      board cpy = *this;
      cpy.white.enpassant = cpy.black.enpassant = 0;
      cpy.white.lastdouble = cpy.black.lastdouble = 0;
      if(m.type == KING) cpy.white.OO = cpy.white.OOO = false;
      if(m.type == ROOK && bool(from & context_OO(white_pov))) cpy.white.OO = false;
      if(m.type == ROOK && bool(from & context_OOO(white_pov))) cpy.white.OOO = false;
      
      cpy.white[m.type] &= ~from;
      cpy.white.all &= ~from;
      cpy.white[m.type] |= to;
      cpy.white.all |= to;
      cpy.black.apply([&to](U64& b){
        b &= ~to;
      });
      
      cpy.black.OOO &= bool(cpy.black.r & context_OOO(!white_pov));
      cpy.black.OO &= bool(cpy.black.r & context_OO(!white_pov));
      ++cpy.count;     
      return cpy;
    }
    
    board do_(const special_move& in){
      board cpy = *this;
      cpy.white = in.result;
      cpy.black.apply([this, in](U64& channel){ channel &= ~(in.result.all | in.eliminated); });
      cpy.black.enpassant = cpy.black.lastdouble = 0;
      ++cpy.count;
			return cpy;
    }

		board& flip(){
			white_pov = !white_pov;
			swap(white, black);
			return *this;
		}

    board rotate() const {
      auto cpy = *this;
      swap(cpy.white, cpy.black);
      cpy.white.apply(util::rotate<U64>);
      cpy.black.apply(util::rotate<U64>);
      return cpy;
    }

		board(){ initmagicmoves(); }
    
  };
  
  ostream& operator<<(ostream& o, board b){
    const state& black = b.white_pov ? b.black : b.white;
    const state& white = b.white_pov ? b.white : b.black;
    for(size_t i(0); i < 8; ++i){
      for(size_t j(0); j < 8; ++j){
        if(white.k & (U64(1) << (i * 8 + j))){ o << " ♔"; continue; }
        if(white.q & (U64(1) << (i * 8 + j))){ o << " ♕"; continue; }
        if(white.r & (U64(1) << (i * 8 + j))){ o << " ♖"; continue; }
        if(white.b & (U64(1) << (i * 8 + j))){ o << " ♗"; continue; }
        if(white.n & (U64(1) << (i * 8 + j))){ o << " ♘"; continue; }
        if(white.p & (U64(1) << (i * 8 + j))){ o << " ♙"; continue; }
        if(black.k & (U64(1) << (i * 8 + j))){ o << " ♚"; continue; }
        if(black.q & (U64(1) << (i * 8 + j))){ o << " ♛"; continue; }
        if(black.r & (U64(1) << (i * 8 + j))){ o << " ♜"; continue; }
        if(black.b & (U64(1) << (i * 8 + j))){ o << " ♝"; continue; }
        if(black.n & (U64(1) << (i * 8 + j))){ o << " ♞"; continue; }
        if(black.p & (U64(1) << (i * 8 + j))){ o << " ♟"; continue; }
        o << " .";        
      }
      o << '\n';
    }
    return o;
  }
  
}

size_t check_mate_counter = 0;

template<typename T> auto negamax(T&& eval, chess::board b, const size_t& depth) -> decltype(eval(b)){
  if(depth == 0) return b.white_pov ? eval(b) : -eval(b);
  optional<decltype(eval(b))> result{};
  b.visit_positions([&result, &eval, &depth](chess::board b_prime){
    if(!result.has_value()) result = -negamax(forward<T>(eval), b_prime, (depth - 1));
    else result = max(result.value(), -negamax(forward<T>(eval), b_prime, (depth - 1))); 
  });
  if(result.has_value()) return result.value();
  else if(b.is_check()) return -numeric_limits<decltype(eval(b))>::max();
  else return decltype(eval(b))(0);
}

template<typename T> auto flatten(chess::board bd){
  using V = typename T::value_type;
  if(!bd.white_pov) bd.flip();
  array<V, 25> result;
  auto iter = result.begin();
  auto write = [&result, &iter](auto functor){
    size_t i(0);
    for(chess::piece_t p : chess::pieces){
      if(iter == result.end()) break;
      *iter = V(functor(p));
      ++iter;
    }
  };

  auto pos = [&bd](chess::piece_t i){
    return bd.white[i];
  };
  
  auto att = [&bd](chess::piece_t i){
    chess::U64 u{0};
    for(size_t val : util::serialize(bd.white[i])) u |= bd.standard_moves(i, val);
    return u;
  };
  
  write(pos);
  write(att);
  bd.flip();
  write(pos);
  write(att);
  *iter = V::filled(bd.count);
  return result;
}

template<typename T, size_t depth, bool P = true> int compare(T a, T b){
  
  auto a_eval = [&a](const chess::board& bd){
    auto k = bd.rotate();
    return a(flatten<T>(bd)).value() - a(flatten<T>(k)).value();
  };
  
  auto b_eval = [&b](const chess::board& bd){
    auto k = bd.rotate();
    return b(flatten<T>(bd)).value() - b(flatten<T>(k)).value();
  };
  
  chess::board bd;
  for(size_t i(0); i < 256; ++i){
    vector<chess::board> options;
    bd.visit_positions([&options](auto in){
      options.push_back(in);
    });
    
    if(options.size() == 0){
      if(bd.is_check()){
        if constexpr (P){
          cout << "loss for " << (bd.white_pov ? "a" : "b") << '\n';
          cout << bd << '\n';
        }
        ++check_mate_counter;
        return !bd.white_pov ? 1 : -1;
      }
      else break;
    }
    
    auto [iter, eval] = util::best_element(options.begin(), options.end(),
        [&](chess::board game){
          if(bd.white_pov) return -negamax(a_eval, game, depth);
          else return -negamax(b_eval, game, depth); 
        }); 

    if constexpr(P){
      cout << "eval :: " << eval << '\n';
      cout << bd << '\n';
    }
    
    bd = *iter;
    
    /*if(bd.white_pov){
      if((bd.white == prev.white) && (bd.black == prev.black)) break;
      else prev = bd;
    }*/  
  }
  return 0;
}



template<typename I> struct integral_board{
  static const I max = 1000;
  static const I min = -max;
  using U64 = unsigned long long;
  using value_type = I;
  array<I, 64> bd;
  
  string to_string(){
    string o;
    for(size_t i(0); i < 8; ++i){
      for(size_t j(0); j < 8; ++j){
        o += std::to_string(bd[i * 8 + j]) + " ";
      }
      o += '\n';
    }
    return o;  
  }  
  
  template<typename P> integral_board<I>& permute(P pr){
    for(auto& val : bd){
      if(util::rand_real<P>() < pr) val = (rand() % (max - min)) + min;
    }
    return *this;
  }
  
  double value(){
    double val = 0;
    for(auto& elem : bd) val += elem;
    val /= 64.0;
    return val;
  }
  
  template<typename Input> integral_board<I>& load(Input& in){
    for(I& val : bd) in >> val;
    return *this;
  }
  
  static integral_board<I> filled(I i){ integral_board ret; ret.bd.fill(i); return ret; }
  
  integral_board(){ bd.fill(I(0)); }
  
  integral_board(U64 i){
    for(size_t val(0); val < 64; ++val){
      bd[val] = bool(i & (U64(1) << val)) ? max : 0;
    }
  }
};

template<typename T, typename Sq> T shift(const T& state, Sq s){
  T result;
  for(size_t i(0); i < 64; ++i){
    chess::sq_t inflated = chess::inflate(i).add(s);
    
    if(chess::s_is_valid(inflated)){
      result.bd[chess::deflate(inflated)] = state.bd[i];
    }
  }
  return result;
}

template<typename A, typename B> A convolve(A strip, B kernel){
  //terrible
  A result; result.fill(0);
  for(int i(0); i < 8; ++i){
    for(int j(0); j < 8; ++j){
      for(int k(-1); k < 2 && (j + k) > 0 && (j + k) < 8; ++k){
        for(int l(-1); l < 2 && (i + l) > 0 && (i + l) < 8; ++l){
          result[i * 8 + j] += strip[(i + l) * 8 + (j + k)] * kernel[k + 1][l + 1];
        }
      }
    }
  }
  return result;
}

int main(){
  srand(time(0));
  using B = integral_board<int>;
  using CP = cartes::program<6144, 25, 64, B>;
  double mutation_rate; cout << "mutation_rate :: "; cin >> mutation_rate;
  double replacement_rate; cout << "replacement_rate :: "; cin >> replacement_rate;
  vector<function<B(B, B)>> grammar = {
    [](B a, B b){
      size_t o = distance(b.bd.begin(), max_element(b.bd.begin(), b.bd.end()));
      return shift(a, chess::inflate(o));
    },
    [](B a, B b){
      for(size_t i(0); i < 64; ++i) a.bd[i] = a.bd[i] + b.bd[i];
      return shift(a, chess::sq_t{0, 1});
    },
    [](B a, B b){
      for(size_t i(0); i < 64; ++i) a.bd[i] = a.bd[i] + b.bd[i];
      return shift(a, chess::sq_t{1, 0});
    },
    [](B a, B b){
      for(size_t i(0); i < 64; ++i) a.bd[i] = a.bd[i] + b.bd[i];
      return shift(a, chess::sq_t{0, -1});
    },
    [](B a, B b){
      for(size_t i(0); i < 64; ++i) a.bd[i] = a.bd[i] + b.bd[i];
      return shift(a, chess::sq_t{-1, 0});
    },
    [](B a, B b){
      for(size_t i(0); i < 64; ++i) a.bd[i] = a.bd[i] + b.bd[i];
      reverse(a.bd.begin(), a.bd.end());
      return a;
    },
    [](B a, B b){
      for(size_t i(0); i < 64; ++i) a.bd[i] = a.bd[i] * b.bd[i];
      return a;
    },
    [](B a, B b){
      for(size_t i(0); i < 64; ++i) a.bd[i] = a.bd[i] + b.bd[i];
      return a;
    },
    [](B a, B b){
      for(size_t i(0); i < 64; ++i) a.bd[i] = a.bd[i] - b.bd[i];
      return a;
    },
    [](B a, B b){
      for(size_t i(0); i < 64; ++i) a.bd[i] = (b.bd[i] != 0) ? (a.bd[i] / b.bd[i]) : 0;
      return a;
    },
    [](B a, B b){
      // maybe ???
      for(size_t i(0); i < 64; ++i) a.bd[i] = a.bd[i] + b.bd[i]  + (rand() % 12);
      return a;      
    },
    [](B a, B b){
      for(size_t i(0); i < 64; ++i){
        auto temp = b.bd[i] + a.bd[i];
        if(temp > 0) a.bd[i] = temp;
      }
      return a;
    },
    [](B a, B b){
      for(size_t i(0); i < 64; ++i) a.bd[i] = b.bd[i] + a.bd[i];
      transform(a.bd.begin(), a.bd.end(), a.bd.begin(), [](auto in){
        return (in > 0) ? in % B::max : in;
      });
      return a;
    },
    [](B a, B b){
      transform(a.bd.begin(), a.bd.end(), a.bd.begin(), [](auto in){
        return -in;
      });
      return a;
    },
    [](B a, B b){
      transform(a.bd.begin(), a.bd.end(), a.bd.begin(), [](auto in){
        return abs(in);
      });
      return a;
    },
    [](B a, B b){
      for(size_t i(0); i < 64; ++i) a.bd[i] = (a.bd[i] > 0) ? b.bd[i] : 0;
      return a;
    },
    [](B a, B b){
      for(size_t i(0); i < 64; ++i) a.bd[i] = (a.bd[i] > 0) ? a.bd[i] : 0;
      return a;
    },
    [](B a, B b){
      auto difference = [](auto b1, auto b2){
        int val = 0;
        for(size_t i(0); i < 64; ++i) val += pow(b1.bd[i] - b2.bd[i], 2);
        return val;
      };
      B c;
      for(size_t i(0); i < 64; ++i) c.bd[i] = difference(shift(a, chess::inflate(i)), b);
      return c;
    },
    [](B a, B b){
      auto difference = [](auto b1, auto b2){
        B::value_type val = 0;
        for(size_t i(0); i < 64; ++i) val += pow(b1.bd[i] - b2.bd[i], 2);
        return val;
      };
      int best_val = difference(a, b);
      int best_i = 0;
      for(int i(0); i < 64; ++i){
        auto val = difference(shift(a, chess::inflate(i)), b);
        if(val < best_val){ best_val = val; best_i = i; }
      }
      a.bd.fill(best_val);
      return a;
    },
    [](B a, B b){
      B::value_type val(0);
      for(size_t i(0); i < 64; ++i){
        val += a.bd[i] + b.bd[i];
      }
      val /= 128;
      a.bd.fill(val);
      return a;
    },
    [](B a, B b){
      auto val = [](B in){
        B::value_type result(0);
        for(auto& v : in.bd) result += v;
        return result;
      };
      return (val(a) > val(b)) ? a : b;
    },
    [](B a, B b){
      auto val = [](B in){
        B::value_type result(0);
        for(auto& v : in.bd) result += v;
        return result;
      };
      return (val(a) < val(b)) ? a : b;
    },
    [](B a, B b){
      B::value_type val(0);
      for(auto& v : a.bd){
        if(v > 0) ++val;
      };
      a.bd.fill(val);
      return a;
    },
    [](B a, B b){
      auto val = [](B in){
        B::value_type result(0);
        for(auto& v : in.bd) result += v;
        return result;
      };
      B::value_type temp = val(b);
      for(auto& v : a.bd) v *= temp;
      return a;
    },
    [](B a, B b){
      for(size_t i(0); i < 64; ++i) if(b.bd[i] > 0) a.bd[i] += 1;
      return a;      
    },
    [](B a, B b){
      for(size_t i(0); i < 64; ++i) if(b.bd[i] > 0) a.bd[i] -= 1;
      return a;      
    },
    [](B a, B b){
      return a;
    },
    [](B a, B b){
      return b;
    },
    [](B a, B b){
      auto dot_prod = [](auto _1, auto _2){
        B::value_type result(0);
        for(size_t val(0); val < 64; ++val) result += _1.bd[val] * _2.bd[val];
        return result;
      };
      B c;
      for(size_t val(0); val < 64; ++val) c.bd[val] = dot_prod(shift(a, chess::inflate(val)), b);
      return c;
    },
    [](B a, B b){
      B::value_type avg(0);
      for(auto val : b.bd) avg += val;
      avg /= 64;
      for(auto& val : a.bd) val = (val >= avg) ? avg : B::value_type(0);
      return a;
    },
    [](B a, B b){
      array<array<B::value_type, 3>, 3> edge
      {{
        array<B::value_type, 3>{{-1, -1, -1}},
        array<B::value_type, 3>{{-1, 8, -1}},
        array<B::value_type, 3>{{-1, -1, -1}}
      }};
      b.bd = convolve(a.bd, edge);
      return b;
    },
    [](B a, B b){
      array<array<B::value_type, 3>, 3> blur
      {{
        array<B::value_type, 3>{{1, 2, 1}},
        array<B::value_type, 3>{{2, 4, 2}},
        array<B::value_type, 3>{{1, 2, 1}}
      }};
      b.bd = convolve(a.bd, blur);
      for(auto& val : b.bd) val /= 16;
      return b;
    },
    [](B a, B b){
      array<array<B::value_type, 3>, 3> sharp
      {{
        array<B::value_type, 3>{{0, -1, 0}},
        array<B::value_type, 3>{{-1, 5, -1}},
        array<B::value_type, 3>{{0, -1, 0}}
      }};
      b.bd = convolve(a.bd, sharp);
      return b;
    },
    [](B a, B b){
      array<array<B::value_type, 3>, 3> edge;
      for(size_t i(0); i < 9; ++i) edge[i / 3][i % 3] = b.bd[i];
      b.bd = convolve(a.bd, edge);
      return b;
    },
    [](B a, B b){
      B c;
      for(auto [a_it, b_it] = make_pair(a.bd.begin(), b.bd.begin()); a_it != a.bd.end(); ++a_it, ++b_it){
        c.bd[abs(*a_it) % 64] = *b_it;
      }
      return c;
    },
    [](B a, B b){
      B c;
      for(auto [a_it, b_it] = make_pair(a.bd.begin(), b.bd.begin()); a_it != a.bd.end(); ++a_it, ++b_it){
        *b_it = max(*b_it, *a_it);
      }
      return c;      
    },
    [](B a, B b){
      B c;
      for(auto [a_it, b_it] = make_pair(a.bd.begin(), b.bd.begin()); a_it != a.bd.end(); ++a_it, ++b_it){
        *b_it = min(*b_it, *a_it);
      }
      return c;      
    },
    [](B a, B b){
      sort(b.bd.begin(), b.bd.end());
      return b;
    }
  };
  ga::population<CP, 6> pop(grammar);
  cout << "load from input file (y/n) :: ";
  string reply{};
  if(cin >> reply; reply == "y"){
    cout << "input file :: ";
    string fname{}; cin >> fname;
    fstream input_file(fname);
    pop.seed = CP::load(grammar, input_file);
    cout << pop.seed.to_string() << '\n';
    pop.evolve(mutation_rate);
  }
  cout << "save file :: ";
  string sf{}; cin >> sf;
  for(size_t i (0); true; ++i){
    if(i % 10 == 0) pop.challenge_seed(compare<CP, 1, true>, 20);
    else pop.challenge_seed(compare<CP, 1, false>, 20);
    fstream(sf, ios_base::out) << pop.seed.to_string() << flush;
    cout << "CHECK_MATE_COUNTER :: " << check_mate_counter << '\n';
    pop.evolve(mutation_rate, replacement_rate);
  }
}
