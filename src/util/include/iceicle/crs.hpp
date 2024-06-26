/**
 * @brief compressed row storage 
 * @author Gianni Absillis (gabsill@ncsu.edu)
 */

#pragma once 

#include <type_traits>
#include <vector>
#include <span>
#include <algorithm>
namespace iceicle::util {

    /**
     * @brief compressed row storage 
     * @tparam T the data type 
     * @tparam IDX the index type 
     * TODO: Allocator
     */
    template<class T, class IDX = std::size_t>
    class crs {
        public:

        // ============
        // = Typedefs =
        // ============
        using value_type = T;
        using index_type = IDX;
        using size_type = std::make_unsigned_t<index_type>;

        private:

        /// @brief the number of nonzero entries (the data ptr size)
        size_type _nnz;

        /// @brief the number of rows
        size_type _nrow;

        /// @brief the data (size = nnz)
        value_type *_data{nullptr};

        /// @brief the indices of the start of each row (size = nrow + 1)
        index_type *_cols{nullptr};

        public:

        // ================
        // = Constructors =
        // ================
        
        constexpr crs() noexcept = default;

        /**
         * @brief construct from existing ragged data 
         * @param ragged_data a 2D ragged array of data to copy
         */
        constexpr crs(std::vector<std::vector<T>> &ragged_data)
        : _nnz{0}, _nrow{ragged_data.size()}, _cols{new index_type[_nrow + 1]}{
            // count up the number of nonzeros and row lengths
            _cols[0] = 0;
            for(index_type irow = 0; irow < ragged_data.size(); ++irow){
                _nnz += ragged_data[irow].size();
                _cols[irow + 1] = _cols[irow] + ragged_data[irow].size();
            }
            _data = new value_type[_nnz];
            for(index_type irow = 0; irow < ragged_data.size(); ++irow){
                std::copy(ragged_data[irow].begin(),
                        ragged_data[irow].end(), _data + _cols[irow]);
            }
        }

        /// @brief copy constructor
        constexpr crs(const crs& other)
        : _nnz(other._nnz), _nrow(other._nrow),
        _data{new value_type[_nnz]}, _cols{new index_type[_nrow + 1]}
        {
            std::copy_n(other._data, _nnz, _data);
            std::copy_n(other._cols, _nrow + 1, _cols);
        }

        /// @brief move constructor
        constexpr crs(crs&& other)
        : _nnz(other._nnz), _nrow(other._nrow), 
        _data{other._data}, _cols{other._cols}
        {
            other._data = nullptr;
            other._cols = nullptr;
        }

        /// @brief copy assignment
        constexpr auto operator=(const crs& other) -> crs<T, IDX>& {
            _nnz = other._nnz;
            _nrow = other._nrow;

            if(_data != nullptr){
                delete[] _data;
            }
            if(_cols != nullptr){
                delete[] _cols;
            }
            std::copy_n(other._data, _nnz, _data);
            std::copy_n(other._cols, _nrow + 1, _cols);
            return *this;
        }

        /// @brief move assignment
        constexpr auto operator=(crs&& other) -> crs<T, IDX>& {
            _nnz = other._nnz;
            _nrow = other._nrow;

            if(_data != nullptr){
                delete[] _data;
            }
            if(_cols != nullptr){
                delete[] _cols;
            }
            _data = other._data;
            other._data = nullptr;
            _cols = other._cols;
            other._cols = nullptr;
            return *this;
        }

        /// @brief destructor
        ~crs(){
            if(_data != nullptr){
                delete[] _data;
            }
            if(_cols != nullptr){
                delete[] _cols;
            }
        }

        // =========
        // = Sizes =
        // =========

        /// @brief get the total number of values stored 
        /// abbreviation for "number of non-zeros" from sparse matrix terminology
        inline constexpr 
        auto nnz() const noexcept -> size_type { return _nnz; }

        //// @brief the number of rows represented
        inline constexpr 
        auto nrow() const noexcept -> size_type { return _nrow; }

        // ============
        // = Indexing =
        // ============
       
        /**
         * @brief index using a 2d index 
         * @param irow the row index
         * @param jcol the column index
         * @return the value at the location of the index pair 
         * (irow, jcol)
         */
        inline constexpr
        auto operator[](index_type irow, index_type jcol) noexcept 
        -> value_type& {
            return _data[_cols[irow] + jcol];
        }
        
        /**
         * @brief index using a 2d index 
         * @param irow the row index
         * @param jcol the column index
         * @return the value at the location of the index pair 
         * (irow, jcol)
         */
        inline constexpr
        auto operator[](index_type irow, index_type jcol) const noexcept 
        -> const value_type& {
            return _data[_cols[irow] + jcol];
        }

        /**
         * @brief get a std::span for a given row 
         * @param irow the index of the row 
         * @return a std::span covering the range of values in the row
         */
        inline constexpr 
        auto rowspan(index_type irow) noexcept 
        -> std::span<value_type> {
            return std::span<value_type>{
                _data + _cols[irow],
                _data + _cols[irow + 1]
            };
        }

        /**
         * @brief get a std::span for a given row 
         * @param irow the index of the row 
         * @return a std::span covering the range of values in the row
         */
        inline constexpr 
        auto rowspan(index_type irow) const noexcept
        -> std::span<const value_type> {
            return std::span<value_type>{
                _data + _cols[irow],
                _data + _cols[irow + 1]
            };
        }
    };

    template<class T>
    crs(std::vector<std::vector<T>>) -> crs<T>;
}
