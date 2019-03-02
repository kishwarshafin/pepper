import torch.nn as nn
import torch.nn.functional as F


def conv3x3(in_planes, out_planes, stride=1):
    """3x3 convolution with padding"""
    return nn.Conv2d(in_planes, out_planes, kernel_size=3, stride=stride,
                     padding=1, bias=False)


class BasicConv2d(nn.Module):

    def __init__(self, in_channels, out_channels, **kwargs):
        super(BasicConv2d, self).__init__()
        self.conv = nn.Conv2d(in_channels, out_channels, bias=False, **kwargs)
        self.bn = nn.BatchNorm2d(out_channels, eps=0.001)

    def forward(self, x):
        x = self.conv(x)
        x = self.bn(x)
        return F.relu(x, inplace=True)


class BasicBlock(nn.Module):
    expansion = 1

    def __init__(self, inplanes, planes, stride=1, downsample=None):
        super(BasicBlock, self).__init__()
        self.conv1 = conv3x3(inplanes, planes, stride)
        self.bn1 = nn.BatchNorm2d(planes)
        self.relu = nn.ReLU(inplace=True)
        self.conv2 = conv3x3(planes, planes)
        self.bn2 = nn.BatchNorm2d(planes)
        self.downsample = downsample
        self.stride = stride

    def forward(self, x):
        residual = x

        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)

        out = self.conv2(out)
        out = self.bn2(out)

        if self.downsample is not None:
            residual = self.downsample(x)

        out += residual
        out = self.relu(out)

        return out


class Bottleneck(nn.Module):
    expansion = 4

    def __init__(self, inplanes, planes, stride=1, downsample=None):
        super(Bottleneck, self).__init__()
        self.conv1 = nn.Conv2d(inplanes, planes, kernel_size=1, bias=False)
        self.bn1 = nn.BatchNorm2d(planes)
        self.conv2 = nn.Conv2d(planes, planes, kernel_size=3, stride=stride,
                               padding=1, bias=False)
        self.bn2 = nn.BatchNorm2d(planes)
        self.conv3 = nn.Conv2d(planes, planes * self.expansion, kernel_size=1, bias=False)
        self.bn3 = nn.BatchNorm2d(planes * self.expansion)
        self.relu = nn.ReLU(inplace=True)
        self.downsample = downsample
        self.stride = stride

    def forward(self, x):
        residual = x

        out = self.conv1(x)
        out = self.bn1(out)
        out = self.relu(out)

        out = self.conv2(out)
        out = self.bn2(out)
        out = self.relu(out)

        out = self.conv3(out)
        out = self.bn3(out)

        if self.downsample is not None:
            residual = self.downsample(x)

        out += residual
        out = self.relu(out)

        return out


class ResNet(nn.Module):

    def __init__(self, in_channels, block, layers):
        self.inplanes = 80
        super(ResNet, self).__init__()
        self.Context_Conv2d_0a = BasicConv2d(in_channels, 14, kernel_size=(1, 3), padding=(0, 1), groups=in_channels)
        self.Context_Conv2d_0b = BasicConv2d(14, 28, kernel_size=(1, 3), padding=(0, 1), stride=(1, 2), groups=14)
        self.Context_Conv2d_0c = BasicConv2d(28, 56, kernel_size=(1, 3), padding=(0, 1), groups=28)
        self.Conv2d_1a_3x3 = BasicConv2d(56, 80, kernel_size=(1, 3), padding=0, stride=(1, 2))

        # self.Context_Conv2d_0a = BasicConv2d(in_channels, 20, kernel_size=3, padding=(1, 0), groups=in_channels)
        # self.Context_Conv2d_0b = BasicConv2d(20, 40, kernel_size=3, padding=(1, 0), groups=20)
        # self.Context_Conv2d_0c = BasicConv2d(40, 80, kernel_size=3, padding=(1, 0), groups=40)
        # self.Conv2d_1a_3x3 = BasicConv2d(80, 80, kernel_size=3, padding=(1, 0))

        self.layer1 = self._make_layer(block, 128, layers[0])
        self.layer2 = self._make_layer(block, 192, layers[1], stride=(1, 2))
        self.layer3 = self._make_layer(block, 256, layers[2], stride=(1, 2))
        self.layer4 = self._make_layer(block, 512, layers[3], stride=(1, 2))

        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
            elif isinstance(m, nn.BatchNorm2d):
                nn.init.constant_(m.weight, 1)
                nn.init.constant_(m.bias, 0)

    def _make_layer(self, block, planes, blocks, stride=(1, 1)):
        downsample = None
        if stride != 1 or self.inplanes != planes * block.expansion:
            downsample = nn.Sequential(
                nn.Conv2d(self.inplanes, planes * block.expansion,
                          kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm2d(planes * block.expansion),
            )
        layers = []
        layers.append(block(self.inplanes, planes, stride, downsample))
        self.inplanes = planes * block.expansion
        for i in range(1, blocks):
            layers.append(block(self.inplanes, planes))

        return nn.Sequential(*layers)

    def forward(self, x):
        x = self.Context_Conv2d_0a(x)
        x = self.Context_Conv2d_0b(x)
        x = self.Context_Conv2d_0c(x)
        x = self.Conv2d_1a_3x3(x)

        x = self.layer1(x)
        x = self.layer2(x)
        x = self.layer3(x)
        x = self.layer4(x)

        return x


def resnet18_custom(input_channels):
    """Constructs a ResNet-18 model.
    """
    model = ResNet(input_channels, BasicBlock, [2, 3, 3, 2])

    return model
