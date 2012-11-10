

namespace functional
{

    template<typename ReturnType, typename... Args>
        std::function<ReturnType (Args...)> memoize(std::function<ReturnType (Args...)> func)
        {
            std::map<std::tuple<Args...>, ReturnType> cache;
            return ([=](Args... args) mutable  {
                    std::tuple<Args...> t(args...);
                    if (cache.find(t) == cache.end())                
                        cache[t] = func(args...);
                    return cache[t];
                    });
        }

}
