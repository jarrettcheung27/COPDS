#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <unordered_map>
#include <memory>
#include <string>
#include <vector>

#include "BCJRQSPA.h"
#include "bch_codec.h"

namespace py = pybind11;

namespace
{
	int info_bits_len(const BCJRQSPA &codec)
	{
		return codec.get_block_length() * codec.m_degree;
	}

	int code_bits_len(const BCJRQSPA &codec)
	{
		return codec.get_total_length() * codec.m_degree;
	}

	void require_1d(const py::buffer_info &info, const char *name)
	{
		if (info.ndim != 1)
		{
			throw std::runtime_error(std::string(name) + " must be 1-D");
		}
	}

	void require_2d(const py::buffer_info &info, const char *name)
	{
		if (info.ndim != 2)
		{
			throw std::runtime_error(std::string(name) + " must be 2-D");
		}
	}

	BinaryBCHCodec &get_bch_codec(int n, int k)
	{
		static std::unordered_map<std::string, std::unique_ptr<BinaryBCHCodec>> cache;
		const std::string key = std::to_string(n) + "_" + std::to_string(k);
		auto it = cache.find(key);
		if (it == cache.end())
		{
			auto codec = std::make_unique<BinaryBCHCodec>(n, k);
			it = cache.emplace(key, std::move(codec)).first;
		}
		return *(it->second);
	}
}

PYBIND11_MODULE(fftqspa, m)
{
	m.doc() = "FFTQSPA BiBo encoder/decoder bindings";

	py::class_<BCJRQSPA>(m, "BCJRQSPA")
		.def(py::init<std::string, int, std::string>(), py::arg("parity_filename"), py::arg("max_iteration"), py::arg("mapping_filename"))
		.def("info_bits_len", [](const BCJRQSPA &self)
			 { return info_bits_len(self); })
		.def("code_bits_len", [](const BCJRQSPA &self)
			 { return code_bits_len(self); })
		.def("encoder4bibo", [](BCJRQSPA &self, py::array_t<int, py::array::c_style | py::array::forcecast> info_bits)
			 {
			auto in_info = info_bits.request();
			require_1d(in_info, "info_bits");

			const int expected = info_bits_len(self);
			if (in_info.size != expected) {
				throw std::runtime_error("info_bits length mismatch: expected " + std::to_string(expected));
			}

			py::array_t<int> code_bits(code_bits_len(self));
			auto out_info = code_bits.request();

			self.encoder4BiBo(static_cast<int *>(in_info.ptr), static_cast<int *>(out_info.ptr));
			return code_bits; }, py::arg("info_bits"), R"doc(
Encode binary info bits to binary codeword bits.

Parameters
----------
info_bits : 1-D int array
    Length must be info_bits_len(). Values should be 0/1.

Returns
-------
code_bits : 1-D int array
    Length is code_bits_len(). Values are 0/1.
)doc")
		.def("decode4bibo", [](BCJRQSPA &self, py::array_t<double, py::array::c_style | py::array::forcecast> rr_bits_prob)
			 {
			auto in_info = rr_bits_prob.request();
			require_1d(in_info, "rr_bits_prob");

			const int expected = code_bits_len(self);
			if (in_info.size != expected) {
				throw std::runtime_error("rr_bits_prob length mismatch: expected " + std::to_string(expected));
			}

			py::array_t<int> decoded_bits(expected);
			auto out_info = decoded_bits.request();

			int iter = self.FFTQSPA4BiBo(static_cast<double *>(in_info.ptr), static_cast<int *>(out_info.ptr));
			return py::make_tuple(decoded_bits, iter); }, py::arg("rr_bits_prob"), R"doc(
Decode using FFTQSPA4BiBo.

Parameters
----------
rr_bits_prob : 1-D float array
    Length must be code_bits_len(). Each value is P(bit=0).

Returns
-------
(decoded_bits, iter) : (1-D int array, int)
    decoded_bits length is code_bits_len().
)doc");

	m.def("bch_encode", [](int n, int k, py::array_t<int, py::array::c_style | py::array::forcecast> msgs)
		  {
		auto info = msgs.request();
		require_2d(info, "msgs");
		const int rows = static_cast<int>(info.shape[0]);
		const int cols = static_cast<int>(info.shape[1]);
		if (cols != k) {
			throw std::runtime_error("msgs second dimension must be k");
		}

		auto &codec = get_bch_codec(n, k);
		py::array_t<int> out({rows, n});
		auto out_info = out.request();

		const int *in_ptr = static_cast<int *>(info.ptr);
		int *out_ptr = static_cast<int *>(out_info.ptr);

		for (int i = 0; i < rows; ++i) {
			std::vector<uint8_t> msg(k, 0);
			for (int j = 0; j < k; ++j) {
				msg[j] = static_cast<uint8_t>(in_ptr[i * cols + j] & 1);
			}
			auto enc = codec.encode(msg);
			for (int j = 0; j < n; ++j) {
				out_ptr[i * n + j] = static_cast<int>(enc[j]);
			}
		}
		return out; }, py::arg("n"), py::arg("k"), py::arg("msgs"));

	m.def("bch_decode", [](int n, int k, py::array_t<int, py::array::c_style | py::array::forcecast> msgs)
		  {
		auto info = msgs.request();
		require_2d(info, "msgs");
		const int rows = static_cast<int>(info.shape[0]);
		const int cols = static_cast<int>(info.shape[1]);
		if (cols != n) {
			throw std::runtime_error("msgs second dimension must be n");
		}

		auto &codec = get_bch_codec(n, k);
		py::array_t<int> out({rows, k + 1});
		auto out_info = out.request();

		const int *in_ptr = static_cast<int *>(info.ptr);
		int *out_ptr = static_cast<int *>(out_info.ptr);

		for (int i = 0; i < rows; ++i) {
			std::vector<uint8_t> rx(n, 0);
			for (int j = 0; j < n; ++j) {
				rx[j] = static_cast<uint8_t>(in_ptr[i * cols + j] & 1);
			}
			auto dec = codec.decode(rx);
			out_ptr[i * (k + 1)] = dec.second;
			for (int j = 0; j < k; ++j) {
				out_ptr[i * (k + 1) + (j + 1)] = static_cast<int>(dec.first[j]);
			}
		}
		return out; }, py::arg("n"), py::arg("k"), py::arg("msgs"));

	m.def("bch_decode_and_vote", [](int n1, int k1, int n2, int k2, int n0, py::array_t<int, py::array::c_style | py::array::forcecast> rx_idx, py::array_t<int, py::array::c_style | py::array::forcecast> rx_data)
		  {
		auto idx_info = rx_idx.request();
		auto data_info = rx_data.request();
		require_2d(idx_info, "rx_idx");
		require_2d(data_info, "rx_data");
		const int rows = static_cast<int>(idx_info.shape[0]);
		if (static_cast<int>(data_info.shape[0]) != rows) {
			throw std::runtime_error("rx_idx and rx_data must have same number of rows");
		}
		if (static_cast<int>(idx_info.shape[1]) != n1) {
			throw std::runtime_error("rx_idx second dimension must be n1");
		}
		if (static_cast<int>(data_info.shape[1]) != n2) {
			throw std::runtime_error("rx_data second dimension must be n2");
		}

		auto &codec_idx = get_bch_codec(n1, k1);
		auto &codec_data = get_bch_codec(n2, k2);

		const int *idx_ptr = static_cast<int *>(idx_info.ptr);
		const int *data_ptr = static_cast<int *>(data_info.ptr);

		std::vector<int> counts(n0, 0);
		std::vector<int> ones(static_cast<size_t>(n0) * k2, 0);

		for (int i = 0; i < rows; ++i) {
			std::vector<uint8_t> rx1(n1, 0);
			std::vector<uint8_t> rx2(n2, 0);
			for (int j = 0; j < n1; ++j) {
				rx1[j] = static_cast<uint8_t>(idx_ptr[i * n1 + j] & 1);
			}
			for (int j = 0; j < n2; ++j) {
				rx2[j] = static_cast<uint8_t>(data_ptr[i * n2 + j] & 1);
			}

			auto dec_id = codec_idx.decode(rx1);
			auto dec_data = codec_data.decode(rx2);

			int idx = 0;
			for (int b = 0; b < k1; ++b) {
				idx = (idx << 1) | (dec_id.first[b] & 1);
			}

			if (idx >= 0 && idx < n0 && static_cast<int>(dec_data.first.size()) == k2) {
				counts[idx] += 1;
				for (int b = 0; b < k2; ++b) {
					ones[static_cast<size_t>(idx) * k2 + b] += (dec_data.first[b] & 1);
				}
			}
		}

		py::array_t<double> v_score({n0, k2});
		auto v_info = v_score.request();
		double *v_ptr = static_cast<double *>(v_info.ptr);

		for (int i = 0; i < n0; ++i) {
			if (counts[i] == 0) {
				for (int b = 0; b < k2; ++b) {
					v_ptr[i * k2 + b] = 0.5 + 1e-2;
				}
			} else {
				for (int b = 0; b < k2; ++b) {
					v_ptr[i * k2 + b] = static_cast<double>(ones[static_cast<size_t>(i) * k2 + b]) / static_cast<double>(counts[i]);
				}
			}
		}

		return v_score; }, py::arg("n1"), py::arg("k1"), py::arg("n2"), py::arg("k2"), py::arg("n0"), py::arg("rx_idx"), py::arg("rx_data"));
}
