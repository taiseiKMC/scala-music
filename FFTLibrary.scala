import java.io.File;
import java.io.StringWriter;
import java.io.PrintWriter;
import java.net.URL;
import java.util.Arrays;
import java.util.ResourceBundle;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.IntStream;

import org.apache.commons.math3.complex.Complex;

object FFTLibrary {
    private def copyComplex(x:Complex):Complex
     = new Complex(x.getReal, x.getImaginary)
    private def FFTWrap(x:Array[Complex], inverse:Boolean) : Array[Complex] =
    {
        val n=x.length
        if(n==1)
            return Array(copyComplex(x(0)))

        val even=Range(0,n/2).map(i=>x(i*2)).toArray
        val ev = FFTWrap(even, inverse)

        val odd=Range(0,n/2).map(i=>x(i*2+1)).toArray
        val od = FFTWrap(odd, inverse)

        val inv=if(inverse) 1 else -1
        
        return Range(0,n).map(
            i => ev(i%(n/2)) add (od(i%(n/2)) multiply new Complex(Math.cos(2*Math.PI*i/n), inv*Math.sin(2*Math.PI*i/n)))
            ).toArray
    }
    private def leastPow2(n:Int):Int =
    {
        var t=1
        while(t<n)
            t<<=1
        return t
    }
    private def checkPow2(ary:Array[Complex]):Array[Complex] = 
    {
        val n=leastPow2(ary.length)
        if(ary.length==n) return ary

        val zero = new Complex(0)
        return Range(0,n).map(i=>if(i<ary.length) ary(i) else zero).toArray
    }

    def FFT(x:Array[Double]):Array[Complex] =
    {
        return FFTWrap(checkPow2(x.map(new Complex(_))), false)
    }
    def FFT(x:Array[Complex]):Array[Complex] =
    {
        return FFTWrap(checkPow2(x), false)
    }
    def IFFT(x:Array[Double]):Array[Complex] =
    {
        var ret=FFTWrap(checkPow2(x.map(new Complex(_))), true)

        for(i<- 0 until ret.length)
            ret(i)=ret(i) divide ret.length
        return ret
    }
    def IFFT(x:Array[Complex]):Array[Complex] =
    {
        var ret=FFTWrap(checkPow2(x), true)
        for(i<-0 until ret.length)
            ret(i)=ret(i) divide ret.length
        return ret
    }
    def FourierExtend(x:Array[Complex]):Array[Complex] =
    {
        var ary=Arrays.copyOf(x, leastPow2(x.length))
        for(i<-x.length until ary.length)
            ary(i)=getConj(x(ary.length-i))
        return ary
    }
    def FourierExtend(x:Array[Double]):Array[Double] =
    {
        var ary=Arrays.copyOf(x, leastPow2(x.length))
        for(i<-x.length until ary.length)
            ary(i)=x(ary.length-i)
        return ary
    }
    def getConj(x:Complex):Complex
     = new Complex(x.getReal(), -x.getImaginary())

    def ISTFT(x:Array[Array[Complex]], window:Array[Double], shiftSize:Int, length:Int):Array[Double] =
    {
        for(i<-0 until x.length)
            x(i) = IFFT(x(i))
        
        var wave = new Array[Double](length)
        var wavea = new Array[Double](length)

        for(i<-0 until x.length; j<-0 until x(i).length if j+shiftSize*i < length)
        {
            if(wavea(j+shiftSize*i) < window(j))
            {
                wave(j+shiftSize*i)=x(i)(j).getReal()/window(j)
                wavea(j+shiftSize*i)=window(j)
            }
        }
        return wave
    }
    def ISTFT(x:Array[Array[Complex]], window:Array[Double], shiftSize:Int):Array[Double]
     = ISTFT(x, window, shiftSize, (x.length-1)*shiftSize+x(0).length)
}
