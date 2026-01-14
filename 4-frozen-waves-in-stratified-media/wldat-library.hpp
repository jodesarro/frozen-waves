#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#include <regex>
#include <algorithm>
#include <limits>

void _string_find_replace(std::string &_the_string, const std::string &_string_to_search, const std::string &_string_to_replace)
{
    size_t pos = 0;
    while ((pos = _the_string.find(_string_to_search, pos)) != std::string::npos)
    {
        _the_string.replace(pos, _string_to_search.length(), _string_to_replace);
        pos += _string_to_replace.length();
    }
}

template<typename T>
std::complex<T> _wldata_to_cpp_complex( std::string& _wldata, T _nothing, bool _flags )
{
    // Define wldata
    std::string wldata = _wldata;

    // Remove spaces
    wldata.erase( std::remove_if( wldata.begin(), wldata.end(), ::isspace ), wldata.end() );

    // To lowercase
    std::transform(wldata.begin(), wldata.end(), wldata.begin(),[](unsigned char c){ return std::tolower(c); });

    // C++ standards
    _string_find_replace(wldata, "*^-", "e-");
    _string_find_replace(wldata, "*^+", "e+");
    _string_find_replace(wldata, "*^", "e+");
    _string_find_replace(wldata, "infinity", "inf");
    _string_find_replace(wldata, "indeterminate", "nan");

    // Regex patterns
    std::regex re_real(R"(([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?|[+-]?inf|[+-]?nan))"); // a, +a, -a
    std::regex re_imag(R"(([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?|[+-]?inf|[+-]?nan)\*i)"); // b*i, +b*i, -b*i
    std::regex re_i(R"(([+-]?)i)"); // i, +i, -i
    std::regex re_a_plus_i(R"(([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?|[+-]?inf|[+-]?nan)\+i)"); // a+i, +a+i, -a+i
    std::regex re_a_minus_i(R"(([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?|[+-]?inf|[+-]?nan)\-i)"); // a-i, +a-i, -a-i
    std::regex re_a_plus_bi(R"(([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?|[+-]?inf|[+-]?nan)\+([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?|[+-]?inf|[+-]?nan)\*i)"); // a+b*i
    std::regex re_a_minus_bi(R"(([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?|[+-]?inf|[+-]?nan)\-([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?|[+-]?inf|[+-]?nan)\*i)"); // a-b*i

    auto to_typenamet = [](const std::string& _string, bool _flags) -> T {
        if ( _string == "inf" || _string == "+inf" || _string == "-inf" )
        {
            if ( _flags )
            {
                std::cerr << "[inf] Flag -> Data contains an item of infinity value." << std::endl;
            }
            return std::numeric_limits<T>::infinity();
        }
        else if ( _string == "nan" || _string == "+nan" || _string == "-nan" )
        {
            if ( _flags )
            {
                std::cerr << "[nan] Flag -> Data contains an item of not-a-number value." << std::endl;
            }
            return std::numeric_limits<T>::quiet_NaN();
        }
        else
        {
            return T( std::stod(_string) );
        }
    };

    // Match against patterns
    std::smatch match;
    if ( std::regex_match(wldata, match, re_a_plus_bi) )
    {
        return {to_typenamet(match[1], _flags), to_typenamet(match[2], _flags)};
    }
    else if ( std::regex_match(wldata, match, re_a_minus_bi) )
    {
        return {to_typenamet(match[1], _flags), -to_typenamet(match[2], _flags)};
    }
    else if ( std::regex_match(wldata, match, re_a_plus_i) )
    {
        return {to_typenamet(match[1], _flags), T(1)};
    }
    else if ( std::regex_match(wldata, match, re_a_minus_i) )
    {
        return {to_typenamet(match[1], _flags), T(-1)};
    }
    else if ( std::regex_match(wldata, match, re_imag) )
    {
        return {0.0, to_typenamet(match[1], _flags)};
    }
    else if ( std::regex_match(wldata, match, re_i) )
    {
        return {0.0, (match[1] == "-" ? T(-1) : T(1))};
    }
    else if ( std::regex_match(wldata, match, re_real) )
    {
        return {to_typenamet(match[1], _flags), T(0)};
    }
    else
    {
        if ( _flags )
        {
            std::cerr << "[inv] Flag -> Data contains an item of invalid format: \"" << _wldata << "\", which was converted to a not-a-number value." << std::endl;
        }
        return {std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()};
    }
}

template<typename T>
void _wldat_import( std::string _wldat_path, T * _real_array, std::complex<T> * _complex_array,
    std::vector<std::vector<std::vector<T>>> &_real_vector, std::vector<std::vector<std::vector<std::complex<T>>>> &_complex_vector,
    bool _is_complex, bool _is_vector, int _imax, int _jmax, int _kmax, bool _flags )
{    
    std::fstream wldat_file;
    wldat_file.open( _wldat_path, std::fstream::in );
    char ch;
    int i = 0; int j = 0; int k = 0;
    std::string the_data = "";
    bool ignore_char = false;
    while ( wldat_file.get(ch) )
    {
        // Dealing with the char.
        if ( ch == '{' || ch == '}' || ch == ' ' || ch == '\n' || ch == '\t' )
        {
            ignore_char = true;
        }
        else if ( ch == '(' )
        {
            // A possible start of a comment.
            wldat_file.get(ch);
            if ( ch == '*' )
            {
                // It is the start of a comment.
                bool the_end_of_the_comments = false;
                int comments = 1;
                // Finding the end of the comment.
                while ( wldat_file.get(ch) && !the_end_of_the_comments )
                {
                    if ( ch == '(' )
                    {
                        // A possible nested comment.
                        wldat_file.get(ch);
                        if ( ch == '*' )
                        {
                            // A nested comment.
                            comments += 1;
                        }
                        //wldat_file.unget();
                        //ch == '(';
                    }
                    else if ( ch == '*' )
                    {
                        // A possible end of a (nested or not) comment.
                        wldat_file.get(ch);
                        if ( ch == ')' )
                        {
                            // End of a (nested or not) comment.
                            comments -= 1 ;
                            if ( comments == 0 )
                            {
                                the_end_of_the_comments = true;
                            }
                        }
                    }
                }
                wldat_file.unget();
                ch == ')';
                ignore_char = true;
            }
            else
            {
                // It is not the start of a comment,
                //  treating '(' as any character.
                wldat_file.unget();
                ch == '(';
            }
        }
        else if ( ch == ',' )
        {
            // Storing the the_data.
            if ( _is_complex )
            {
                if ( _is_vector )
                {
                    _complex_vector.at(i).at(j).at(k) = _wldata_to_cpp_complex( the_data, T(0), _flags );
                }
                else
                {
                    _complex_array[i + _imax*(j + _jmax*k)] = _wldata_to_cpp_complex( the_data, T(0), _flags );
                }
            }
            else
            {
                if ( _is_vector )
                {
                    _real_vector.at(i).at(j).at(k) = std::real( _wldata_to_cpp_complex( the_data, T(0), _flags ) );
                }
                else
                {
                    _real_array[i + _imax*(j + _jmax*k)] = std::real( _wldata_to_cpp_complex( the_data, T(0), _flags ) );
                }
            }
            the_data = "";

            // Iterating the indices
            k++;
            if ( k == _kmax )
            {
                k = 0;
                j++;
                if ( j == _jmax )
                {
                    j = 0;
                    i++;
                    if ( i == _imax )
                    {
                        // It shoud be an error.
                        i = j = k = 0;
                    }
                }
            }
            ignore_char = true;
        }
        
        if ( !ignore_char )
        {
            // Adding char to the_data string.
            the_data += std::string(1, ch);
        }
        else
        {
            ignore_char = false;
        }
    }

    // Storing the last the_data of the file.
    if ( _is_complex )
    {
        if ( _is_vector )
        {
            _complex_vector.at(i).at(j).at(k) = _wldata_to_cpp_complex( the_data, T(0), _flags );
        }
        else
        {
            _complex_array[i + _imax*(j + _jmax*k)] = _wldata_to_cpp_complex( the_data, T(0), _flags );
        }
    }
    else
    {
        if ( _is_vector )
        {
            _real_vector.at(i).at(j).at(k) = std::real( _wldata_to_cpp_complex( the_data, T(0), _flags ) );
        }
        else
        {
            _real_array[i + _imax*(j + _jmax*k)] = std::real( _wldata_to_cpp_complex( the_data, T(0), _flags ) );
        }
    }

    wldat_file.close();
}

template<typename T>
void _wldat_export( std::string _wldat_path, T * _real_array, std::complex<T> * _complex_array,
    std::vector<std::vector<std::vector<T>>> &_real_vector, std::vector<std::vector<std::vector<std::complex<T>>>> &_complex_vector,
    bool _is_complex, bool _is_vector, int _imax, int _jmax, int _kmax, int _out_precision , bool _out_scientific, bool _flags )
{
    std::fstream wldat_file;
    wldat_file.open( _wldat_path, std::fstream::out );
    wldat_file << "(*" << " Created with C/C++ using the WLdat library: https://github.com/jodesarro/wldat-library " << "*)" << std::endl;
    wldat_file << "{";
    for ( int i=0; i<_imax; i++ )
    {
        wldat_file << "{";
        for ( int j=0; j<_jmax; j++ )
        {
            wldat_file << "{";
            for ( int k=0; k<_kmax; k++ )
            {
                std::stringstream the_data;
                if ( _out_precision != 0 )
                {
                    the_data.precision(_out_precision);
                }
                if ( _out_scientific )
                {
                    the_data << std::scientific;
                }
                the_data << std::showpos;
                if ( _is_complex )
                {
                    if ( _is_vector )
                    {
                        if ( _flags )
                        {
                            if (  std::isnan( std::real(_complex_vector.at(i).at(j).at(k)) ) ||
                            std::isnan( std::imag(_complex_vector.at(i).at(j).at(k)) ) )
                            {
                                std::cerr << "[nan] Flag -> Data contains an item of not-a-number value." << std::endl;
                            }
                            if ( std::isinf( std::real(_complex_vector.at(i).at(j).at(k)) ) ||
                            std::isinf( std::imag(_complex_vector.at(i).at(j).at(k)) ) )
                            {
                                std::cerr << "[inf] Flag -> Data contains an item of infinity value." << std::endl;
                            }
                        }
                        the_data << std::real(_complex_vector.at(i).at(j).at(k)) << std::imag(_complex_vector.at(i).at(j).at(k)) << "*I";
                    }
                    else
                    {
                        if ( _flags )
                        {
                            if ( std::isnan( std::real(_complex_array[i + _imax*(j + _jmax*k)]) ) ||
                            std::isnan( std::imag(_complex_array[i + _imax*(j + _jmax*k)]) ) )
                            {
                                std::cerr << "[nan] Flag -> Data contains an item of not-a-number value." << std::endl;
                            }
                            if ( std::isinf( std::real(_complex_array[i + _imax*(j + _jmax*k)]) ) ||
                            std::isinf( std::imag(_complex_array[i + _imax*(j + _jmax*k)]) ) )
                            {
                                std::cerr << "[inf] Flag -> Data contains an item of infinity value." << std::endl;
                            }
                        }
                        the_data << std::real(_complex_array[i + _imax*(j + _jmax*k)]) << std::imag(_complex_array[i + _imax*(j + _jmax*k)]) << "*I";
                    }
                }
                else
                {
                    if ( _is_vector )
                    {
                        if ( _flags )
                        {
                            if ( std::isnan( _real_vector.at(i).at(j).at(k) ) )
                            {
                                std::cerr << "[nan] Flag -> Data contains an item of not-a-number value." << std::endl;
                            }
                            if ( std::isinf( _real_vector.at(i).at(j).at(k) ) )
                            {
                                std::cerr << "[inf] Flag -> Data contains an item of infinity value." << std::endl;
                            }
                        }
                        the_data << _real_vector.at(i).at(j).at(k);
                    }
                    else
                    {
                        if ( _flags )
                        {
                            if ( std::isnan( _real_array[i + _imax*(j + _jmax*k)] ) )
                            {
                                std::cerr << "[nan] Flag -> Data contains an item of not-a-number value." << std::endl;
                            }
                            if ( std::isinf( _real_array[i + _imax*(j + _jmax*k)] ) )
                            {
                                std::cerr << "[inf] Flag -> Data contains an item of infinity value." << std::endl;
                            }
                        }
                        the_data << _real_array[i + _imax*(j + _jmax*k)];
                    }
                }
                std::string the_data_str = the_data.str();
                _string_find_replace(the_data_str, "e", "*^");
                _string_find_replace(the_data_str, "+nan", "+Indeterminate");
                _string_find_replace(the_data_str, "-nan", "-Indeterminate");
                _string_find_replace(the_data_str, "nan", "+Indeterminate");
                _string_find_replace(the_data_str, "+inf", "+Infinity");
                _string_find_replace(the_data_str, "-inf", "-Infinity");
                _string_find_replace(the_data_str, "inf", "+Infinity");
                wldat_file << the_data_str;
                ( k == _kmax-1 ) ? ( wldat_file << "}") : ( wldat_file << ",");
            }
            ( j == _jmax-1 ) ? ( wldat_file << "}") : ( wldat_file << ",") ;
        }
        ( i == _imax-1 ) ? ( wldat_file << "}") : ( wldat_file << ",") ;
    }
    wldat_file.close();
}

template<typename T>
void wldat_import( std::string wldat_path, T * data_array, int imax, int jmax, int kmax, bool flags = false )
{
    std::complex<T> empty_complex_array [1];
    std::vector<std::vector<std::vector<T>>> empty_real_vector(1, std::vector<std::vector<T>>(1, std::vector<T>(1) ) );
    std::vector<std::vector<std::vector<std::complex<T>>>> empty_complex_vector(1, std::vector<std::vector<std::complex<T>>>(1, std::vector<std::complex<T>>(1) ) );
    _wldat_import( wldat_path, data_array, empty_complex_array, empty_real_vector, empty_complex_vector, false, false, imax, jmax, kmax, flags );
}

template<typename T>
void wldat_import( std::string wldat_path, std::complex<T> * data_array, int imax, int jmax, int kmax, bool flags = false )
{
    T empty_real_array [1];
    std::vector<std::vector<std::vector<T>>> empty_real_vector(1, std::vector<std::vector<T>>(1, std::vector<T>(1) ) );
    std::vector<std::vector<std::vector<std::complex<T>>>> empty_complex_vector(1, std::vector<std::vector<std::complex<T>>>(1, std::vector<std::complex<T>>(1) ) );
    _wldat_import( wldat_path, empty_real_array, data_array, empty_real_vector, empty_complex_vector, true, false, imax, jmax, kmax, flags );
}

template<typename T>
void wldat_import( std::string wldat_path, std::vector<std::vector<std::vector<T>>> &data_vector, int imax, int jmax, int kmax, bool flags = false )
{
    T empty_real_array [1];
    std::complex<T> empty_complex_array [1];
    std::vector<std::vector<std::vector<std::complex<T>>>> empty_complex_vector(1, std::vector<std::vector<std::complex<T>>>(1, std::vector<std::complex<T>>(1) ) );
    _wldat_import( wldat_path, empty_real_array, empty_complex_array, data_vector, empty_complex_vector, false, true, imax, jmax, kmax, flags );
}

template<typename T>
void wldat_import( std::string wldat_path, std::vector<std::vector<std::vector<std::complex<T>>>> &data_vector, int imax, int jmax, int kmax, bool flags = false )
{
    T empty_real_array [1];
    std::complex<T> empty_complex_array [1];
    std::vector<std::vector<std::vector<T>>> empty_real_vector(1, std::vector<std::vector<T>>(1, std::vector<T>(1) ) );
    _wldat_import( wldat_path, empty_real_array, empty_complex_array, empty_real_vector, data_vector, true, true, imax, jmax, kmax, flags );
}

template<typename T>
void wldat_export( std::string wldat_path, T * data_array, int imax, int jmax, int kmax, int out_precision = 0, bool out_scientific = false, bool flags = false )
{
    std::complex<T> empty_complex_array [1];
    std::vector<std::vector<std::vector<T>>> empty_real_vector(1, std::vector<std::vector<T>>(1, std::vector<T>(1) ) );
    std::vector<std::vector<std::vector<std::complex<T>>>> empty_complex_vector(1, std::vector<std::vector<std::complex<T>>>(1, std::vector<std::complex<T>>(1) ) );
    _wldat_export( wldat_path, data_array, empty_complex_array, empty_real_vector, empty_complex_vector, false, false, imax, jmax, kmax, out_precision, out_scientific, flags );
}

template<typename T>
void wldat_export( std::string wldat_path, std::complex<T> * data_array, int imax, int jmax, int kmax, int out_precision = 0, bool out_scientific = false, bool flags = false )
{
    T empty_real_array [1];
    std::vector<std::vector<std::vector<T>>> empty_real_vector(1, std::vector<std::vector<T>>(1, std::vector<T>(1) ) );
    std::vector<std::vector<std::vector<std::complex<T>>>> empty_complex_vector(1, std::vector<std::vector<std::complex<T>>>(1, std::vector<std::complex<T>>(1) ) );
    _wldat_export( wldat_path, empty_real_array, data_array, empty_real_vector, empty_complex_vector, true, false, imax, jmax, kmax, out_precision, out_scientific, flags );
}

template<typename T>
void wldat_export( std::string wldat_path, std::vector<std::vector<std::vector<T>>> &data_vector, int imax, int jmax, int kmax, int out_precision = 0, bool out_scientific = false, bool flags = false )
{
    T empty_real_array [1];
    std::complex<T> empty_complex_array [1];
    std::vector<std::vector<std::vector<std::complex<T>>>> empty_complex_vector(1, std::vector<std::vector<std::complex<T>>>(1, std::vector<std::complex<T>>(1) ) );
    _wldat_export( wldat_path, empty_real_array, empty_complex_array, data_vector, empty_complex_vector, false, true, imax, jmax, kmax, out_precision, out_scientific, flags );
}

template<typename T>
void wldat_export( std::string wldat_path, std::vector<std::vector<std::vector<std::complex<T>>>> &data_vector, int imax, int jmax, int kmax, int out_precision = 0, bool out_scientific = false, bool flags = false )
{
    T empty_real_array [1];
    std::complex<T> empty_complex_array [1];
    std::vector<std::vector<std::vector<T>>> empty_real_vector(1, std::vector<std::vector<T>>(1, std::vector<T>(1) ) );
    _wldat_export( wldat_path, empty_real_array, empty_complex_array, empty_real_vector, data_vector, true, true, imax, jmax, kmax, out_precision, out_scientific, flags );
}

void wldat_getsize( std::string wldat_path, int &imax, int &jmax, int &kmax, bool flags = false )
{
    std::fstream wldat_file;
    wldat_file.open( wldat_path, std::fstream::in );
    char ch;
    char last_ch = '\0';
    int comma = 0;
    int brace1 = 0;
    int brace2 = 0;
    while ( wldat_file >> ch )
    {
        if ( ch == ',' )
        {
            comma++;
        }
        else if ( ch == '{' )
        {
            brace1++;
            if ( last_ch == '{' )
            {
                brace2++;
            }
        }
        last_ch = ch;
    }
    wldat_file.close();
    imax = brace2 - 1;
    jmax = int( (brace1 -  brace2)/imax );
    kmax = int( (comma + 1)/(imax*jmax) );

    if ( flags )
    {
        double empty_real_array [1];
        std::complex<double> * complex_array = new std::complex<double> [imax*jmax*kmax];
        std::vector<std::vector<std::vector<double>>> empty_real_vector(1, std::vector<std::vector<double>>(1, std::vector<double>(1) ) );
        std::vector<std::vector<std::vector<std::complex<double>>>> empty_complex_vector(1, std::vector<std::vector<std::complex<double>>>(1, std::vector<std::complex<double>>(1) ) );
        _wldat_import( wldat_path, empty_real_array, complex_array, empty_real_vector, empty_complex_vector, true, false, imax, jmax, kmax, flags );
    }
}