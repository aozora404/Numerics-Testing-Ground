using ILGPU;
using ILGPU.Runtime;
using System;

static void AddVectorsKernel(Index1 index, ArrayView<float> a, ArrayView<float> b, ArrayView<float> c)
{
    c[index] = a[index] + b[index];
}

static void Main(string[] args)
{
    int length = 1024;
    float[] a = new float[length];
    float[] b = new float[length];
    float[] c = new float[length];

    for (int i = 0; i < length; i++)
    {
        a[i] = i;
        b[i] = 2 * i;
    }

    using var context = new Context();
    using var accelerator = Accelerator.Create(context, context.Accelerators.First());
    using var bufferA = accelerator.Allocate<float>(length);
    using var bufferB = accelerator.Allocate<float>(length);
    using var bufferC = accelerator.Allocate<float>(length);

    bufferA.CopyFrom(a, 0, Index1.Zero, length);
    bufferB.CopyFrom(b, 0, Index1.Zero, length);

    var kernel = accelerator.LoadStreamKernel<ArrayView<float>, ArrayView<float>, ArrayView<float>>(AddVectorsKernel);
    kernel((length, accelerator.DefaultStream), bufferA.View, bufferB.View, bufferC.View);

    accelerator.Synchronize();

    bufferC.CopyTo(c, 0, Index1.Zero, length);

    // Print first 10 results for verification
    for (int i = 0; i < 10; i++)
    {
        Console.WriteLine(c[i]);
    }
}
