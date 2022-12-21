## Rozwiązanie numeryczne równania różniczkowego przy wykorzystaniu MES | Bartosz Hanc | AGH UST 2022/23
### Opis programu
Repozytorium zawiera plik PDF zawierający sformułowanie danego problemu
brzegowego wraz z przekształceniem go do postaci słabej oraz jego dyskretyzację
wykorzystaną w algorytmie numerycznym napisanym w języku Julia. Algorytm
numeryczny wykorzystuje całkowanie numeryczne metodą kwadratur Gaussa-Legendre'a
oraz domyślnie zaimplementowane w Julii, wydajne operacje na macierzach.

### Uruchomienie programu
Do uruchomienia programu wymagany jest zainstalowany [kompilator języka
Julia](https://julialang.org/downloads/) oraz pakiet Plots, który można
zainstalować przez Julia REPL:
```
$ julia
julia> using Pkg
julia> Pkg.add("Plots")
```
Aby uruchomić program, można wykonać polecenie:
```
$ julia solve.jl
```
jednak ze względu na architekturę Julii i długie czasy kompilacji ten sposób nie
jest zalecany. Zalecane jest uruchomienie programu w Julia REPL tak jak w
przykładzie (w przykładzie wprowadzono najpierw liczbę elementów równą 0, aby
skompilować program, a następnie równą 100)
```
$ julia
julia> include("solve.jl")
>> Input number of elements: 0
>> Press enter to exit...

julia> include("solve.jl")
>> Input number of elements: 100
```
