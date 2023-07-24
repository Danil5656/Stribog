#include "GOST3411.h"

namespace Util
{
	//метод для переворачивания последовательности чтобы вывод был корректный
	void reverseVector(uint8_t* vector, uint8_t* vectorReverse, const uint16_t& countElements)
	{
		for (int i = 0; i < countElements; ++i)
		{
			vectorReverse[i] = vector[countElements - 1 - i];
		}
	}

	//метод для вывода полученного хеша на экран
	void printHash(const uint8_t* hash, const uint16_t& sizeHash)
	{
		std::cout << "Hash" << " (" << std::dec << sizeHash << "bits): ";
		for (int i = 0; i < sizeHash / 8; ++i)
		{
			std::cout << std::hex << std::setfill('0') << std::setw(2) << (uint16_t)hash[i];
		}
		std::cout << std::endl;
	}
}

//конструктор, где происходит инициализация полей класса
GOST3411::GOST3411(const uint16_t& sizeHash)
{
	this->sizeHash = sizeHash;
	initializationStage();
}

//деструктор, где происходит освобождения памяти для полей класса
GOST3411::~GOST3411()
{
	delete[] h, checksum, N, bufferDataBlock, hash;
}

//метод, выполняющий побитовое исключающее ИЛИ над 512-битными блоками
void GOST3411::conversionAdd512(uint8_t* a, uint8_t* b, uint8_t* result)
{
	int tmp = 0;
	for (int i = 0; i < this->SIZE_BLOCK_BYTES; i++)
	{
		tmp = a[i] + b[i] + (tmp >> 8);
		result[i] = tmp & 0xff;
	}
}

//метод, выполняющий побитовое сложение по mod 2 над векторами
void GOST3411::conversionX(uint8_t* inputVector, uint8_t* gamma, uint8_t* result)
{
	for (int i = 0; i < this->SIZE_BLOCK_BYTES; ++i)
	{
		result[i] = inputVector[i] ^ gamma[i];
	}
}

//метод, реализующий нелинейное биективное преобразование (подстановку или S-преобразование)
void  GOST3411::conversionS(uint8_t* inputVector)
{
	uint8_t* tmp = new uint8_t[this->SIZE_BLOCK_BYTES]{ 0 };

	for (int i = this->SIZE_BLOCK_BYTES - 1; i >= 0; --i)
	{
		tmp[i] = sboxPi[inputVector[i]];
	}

	std::copy_n(tmp, this->SIZE_BLOCK_BYTES, inputVector);
	delete[] tmp;
}

//метод, реализующий перестановку байт (P-преобразование)
void GOST3411::conversionP(uint8_t* inputVector)
{
	uint8_t* tmp = new uint8_t[this->SIZE_BLOCK_BYTES]{ 0 };

	for (int i = this->SIZE_BLOCK_BYTES - 1; i >= 0; i--)
	{
		tmp[i] = inputVector[tau[i]];
	}

	std::copy_n(tmp, this->SIZE_BLOCK_BYTES, inputVector);
	delete[] tmp;
}

//метод, реализующий линейное преобразование (L-преобразование)
void GOST3411::conversionL(uint8_t* inputVector)
{
	uint64_t* tmpInput = reinterpret_cast<uint64_t*>(inputVector);

	uint8_t indexByte = this->SIZE_BLOCK_BYTES - 1;
	for (int i = 7; i >= 0; --i)
	{
		uint64_t tmpCalcXor = 0;
		for (int j = 63; j >= 0; --j)
		{
			if ((tmpInput[i] >> j) & 1)
			{
				tmpCalcXor ^= matrixA[63 - j];
			}
		}

		for (int j = 0, k = 7; j < 8; j++, --k)
		{
			inputVector[indexByte] = (tmpCalcXor >> k * 8) & 0xff;
			--indexByte;
		}
	}
}

//метод, который получает раундовый ключ и формирует новый ключ
void GOST3411::getRoundKey(uint8_t* K, const int& i)
{
	conversionX(K, const_cast<uint8_t*>(C[i]), K);
	conversionS(K);
	conversionP(K);
	conversionL(K);
}

//метод, реализующий E-преобразование
void GOST3411::conversionE(uint8_t* m, uint8_t* K, uint8_t* result)
{
	conversionX(m, K, result);
	for (int i = 0; i < 12; ++i)
	{
		conversionS(result);
		conversionP(result);
		conversionL(result);
		getRoundKey(K, i);
		conversionX(result, K, result);
	}
}

//метод, реализующий главную функцию сжатия g
void GOST3411::conversionG(uint8_t* h, uint8_t* N, uint8_t* m)
{
	uint8_t* K = new uint8_t[this->SIZE_BLOCK_BYTES]{ 0 };
	uint8_t* tmp = new uint8_t[this->SIZE_BLOCK_BYTES]{ 0 };

	conversionX(h, N, K);
	conversionS(K);
	conversionP(K);
	conversionL(K);

	conversionE(m, K, tmp);
	conversionX(tmp, h, tmp);
	conversionX(tmp, m, h);
	delete[] K, tmp;
}

//метод, реализующий усечение строк длины 512 бит до строки длины 256 бит
void GOST3411::conversionMSB()
{
	for (int i = 0, j = 63; i < this->sizeHash / 8; ++i, --j)
	{
		hash[i] = this->h[j];
	}
}

//метод, реализующий дополнения блока в случае если его длина < 64 байт
void GOST3411::paddingBlock(const uint16_t& sizeBlockBytes, uint8_t* blockData)
{
	if (sizeBlockBytes == this->SIZE_BLOCK_BYTES)
	{
		Util::reverseVector(blockData, this->bufferDataBlock, this->SIZE_BLOCK_BYTES);
		return;
	}
	uint8_t* tmpDataShift = new uint8_t[this->SIZE_BLOCK_BYTES]{ 0 };
	for (int i = this->SIZE_BLOCK_BYTES - 1, j = sizeBlockBytes - 1; j >= 0; --i, --j)
	{
		tmpDataShift[i] = blockData[j];
	}
	tmpDataShift[this->SIZE_BLOCK_BYTES - sizeBlockBytes - 1] = 0x01;
	Util::reverseVector(tmpDataShift, this->bufferDataBlock, this->SIZE_BLOCK_BYTES);
	delete[] tmpDataShift;
}

//метод для инициализации первоначальных значений переменным
void GOST3411::initializationStage()
{
	this->h = new uint8_t[this->SIZE_BLOCK_BYTES]{ 0 };
	this->hash = new uint8_t[this->sizeHash / 8]{ 0 };
	this->N = new uint8_t[this->SIZE_BLOCK_BYTES]{ 0 };
	this->bufferDataBlock = new uint8_t[this->SIZE_BLOCK_BYTES]{ 0 };
	this->checksum = new uint8_t[this->SIZE_BLOCK_BYTES]{ 0 };
	if (this->sizeHash == 512)
	{
		std::copy_n(iv512, this->SIZE_BLOCK_BYTES, this->h);
	}
	else if (this->sizeHash == 256)
	{
		std::copy_n(iv256, this->SIZE_BLOCK_BYTES, this->h);
	}
	iv512[1] = 0x02;
}

//метод по 2 этапу ГОСТ, где хешируем отдельные блоки по 512 бит пока они не кончатся
void GOST3411::secondStage(uint8_t* blockData)
{
	paddingBlock(this->SIZE_BLOCK_BYTES, blockData);
	conversionG(this->h, this->N, this->bufferDataBlock);
	conversionAdd512(this->N, iv512, this->N);
	conversionAdd512(this->checksum, this->bufferDataBlock, this->checksum);
}

//метод по 3 этапу ГОСТ, где хешируем остаток, который не попал в 64-байтный блок
void GOST3411::thirdStage(uint8_t* blockData, const uint16_t& sizeBlockBytes)
{
	uint8_t* tmp = new uint8_t[this->SIZE_BLOCK_BYTES]{ 0 };

	tmp[1] = ((sizeBlockBytes * 8) >> 8) & 0xff;
	tmp[0] = (sizeBlockBytes * 8) & 0xff;

	paddingBlock(sizeBlockBytes, blockData);
	conversionG(this->h, this->N, this->bufferDataBlock);
	conversionAdd512(this->N, tmp, this->N);
	conversionAdd512(this->checksum, this->bufferDataBlock, this->checksum);
	iv512[1] = 0x00;
	conversionG(this->h, this->iv512, this->N);
	conversionG(this->h, this->iv512, this->checksum);
	if (this->sizeHash == 512)
	{
		Util::reverseVector(this->h, this->hash, this->sizeHash / 8);
	}
	else
	{
		conversionMSB();
	}
	Util::printHash(hash, this->sizeHash);
	delete[] tmp;
}

//метод для получения хеш-кода сообщения
void GOST3411::getHash(const uint8_t* blockData, const uint32_t& sizeBlockDataBytes)
{
	if (sizeBlockDataBytes < this->SIZE_BLOCK_BYTES)
	{
		thirdStage(const_cast<uint8_t*>(blockData), sizeBlockDataBytes);
		return;
	}
	else
	{
		uint32_t sizeBlockDataBits = sizeBlockDataBytes * 8;
		uint32_t index = sizeBlockDataBytes - this->SIZE_BLOCK_BYTES;
		while (sizeBlockDataBits >= 512)
		{
			uint8_t* tmpData = new uint8_t[this->SIZE_BLOCK_BYTES]{ 0 };
			for (int i = 0; i < this->SIZE_BLOCK_BYTES; ++i)
			{
				tmpData[i] = blockData[index];
				++index;
			}
			secondStage(tmpData);
			sizeBlockDataBits -= this->SIZE_BLOCK_BYTES * 8;
			index -= this->SIZE_BLOCK_BYTES * 2;
			delete[] tmpData;
		}
		uint8_t* tmpData = new uint8_t[this->SIZE_BLOCK_BYTES]{ 0 };
		for (int i = 0; i < sizeBlockDataBits / 8; ++i)
		{
			tmpData[i] = blockData[i];
		}
		thirdStage(tmpData, sizeBlockDataBits / 8);
		delete[] tmpData;
	}

}